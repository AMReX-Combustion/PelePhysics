target_include_directories(${pelephysics_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(${pelephysics_lib_name}
  PRIVATE
    Constant.H
    Simple.H
    Sutherland.H
    Transport.cpp
    Transport.H
    TransportParams.H
    TransportTypes.H
  )

if("${PELEPHYSICS_TRANSPORT_MODEL}" STREQUAL "Constant")
  target_compile_definitions(${pelephysics_lib_name} PUBLIC USE_CONSTANT_TRANSPORT)
endif()
if("${PELEPHYSICS_TRANSPORT_MODEL}" STREQUAL "Simple")
  target_compile_definitions(${pelephysics_lib_name} PUBLIC USE_SIMPLE_TRANSPORT)
endif()
if("${PELEPHYSICS_TRANSPORT_MODEL}" STREQUAL "Sutherland")
  target_compile_definitions(${pelephysics_lib_name} PUBLIC USE_SUTHERLAND_TRANSPORT)
endif()

