#ifndef TRANSPORTTYPES_H
#define TRANSPORTTYPES_H

namespace pele {
namespace physics {
namespace transport {
struct ConstTransport;
struct SimpleTransport;
struct SutherlandTransport;
} // namespace transport
#ifdef USE_CONSTANT_TRANSPORT
using TransportType = transport::ConstTransport;
#elif USE_SIMPLE_TRANSPORT
using TransportType = transport::SimpleTransport;
#elif USE_SUTHERLAND_TRANSPORT
using TransportType = transport::SutherlandTransport;
#else
static_assert(false, "Invalid Transport specified");
#endif
} // namespace physics
} // namespace pele
#endif
