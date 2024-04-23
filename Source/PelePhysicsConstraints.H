#include <type_traits>

namespace pele::physics {

// Forward declarations
namespace eos {
struct Fuego;
struct GammaLaw;
struct SRK;
struct Manifold;
} // namespace eos

// Forward declarations
namespace transport {
struct ConstTransport;
struct SimpleTransport;
struct SutherlandTransport;
struct ManifoldTransport;
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

// Manifold Transport doesn't apply except with Manifold EOS
// Manifold EOS requires Manifold or Constant Transport
template<typename EosModel>
struct is_valid_physics_combination<EosModel, transport::ManifoldTransport>
  : public std::false_type
{
};

template<typename TransportModel>
struct is_valid_physics_combination<eos::Manifold, TransportModel>
  : public std::false_type
{
};

template<>
struct is_valid_physics_combination<eos::Manifold, transport::ManifoldTransport>
  : public std::true_type
{
};

template<>
struct is_valid_physics_combination<eos::Manifold, transport::ConstTransport>
  : public std::true_type
{
};

} // namespace pele::physics
