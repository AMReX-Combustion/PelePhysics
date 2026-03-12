#ifndef TRANSPORTTYPES_H
#define TRANSPORTTYPES_H

namespace pele::physics {
namespace transport {
struct ConstTransport;
struct SimpleTransport;
struct SutherlandTransport;
struct ManifoldTransport;
} // namespace transport
#ifdef USE_CONSTANT_TRANSPORT
using TransportType = transport::ConstTransport;
#elif USE_SIMPLE_TRANSPORT
using TransportType = transport::SimpleTransport;
#elif USE_SUTHERLAND_TRANSPORT
using TransportType = transport::SutherlandTransport;
#elif USE_MANIFOLD_TRANSPORT
using TransportType = transport::ManifoldTransport;
#else
static_assert(false, "Invalid Transport specified");
#endif
} // namespace pele::physics
#endif
