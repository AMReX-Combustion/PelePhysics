cmake_minimum_required(VERSION 3.23)

project( PelePhysics
    DESCRIPTION "A respository of physics databases and implementation code for use with the Pele suite of of codes"
    HOMEPAGE_URL "https://github.com/AMReX-Combustion/PelePhysics"
    LANGUAGES CXX
    )

message(STATUS "CMake version: ${CMAKE_VERSION}")

#
# Options for model selection
#
set( PELEPHYSICS_SPACEDIM        3        CACHE STRING "Spatial dimension")
set( PELEPHYSICS_PRECISION       "DOUBLE" CACHE STRING "Floating point precision")
set( PELEPHYSICS_CHEMISTRY_MODEL drm19    CACHE STRING "Chemical mechanism")
set( PELEPHYSICS_EOS_MODEL       Fuego    CACHE STRING "Eq. of state model")
set( PELEPHYSICS_TRANSPORT_MODEL Simple   CACHE STRING "Transport model")

#
# Options for library content
#
option( PELEPHYSICS_REACTIONS "Enable Reaction"  ON)
option( PELEPHYSICS_TRANSPORT "Enable Transport" ON)
option( PELEPHYSICS_EOS       "Enable EOS"       ON)
option( PELEPHYSICS_UTILS     "Enable Utilities" ON)

#
# Options for HPC
#
option( PELEPHYSICS_MPI  "Enable MPI"    OFF)
option( PELEPHYSICS_OMP  "Enable OpenMP" OFF)
option( PELEPHYSICS_CUDA "Enable CUDA"   OFF)
option( PELEPHYSICS_HIP  "Enable HIP"    OFF)
option( PELEPHYSICS_SYCL "Enable SYCL"   OFF)
set(PELEPHYSICS_GPU_BACKEND NONE)

if (PELEPHYSICS_CUDA)
   enable_language(CUDA)
   set(PELEPHYSICS_GPU_BACKEND CUDA)
endif ()

if (PELEPHYSICS_HIP)
    enable_language(HIP)
    set(PELEPHYSICS_GPU_BACKEND HIP)
endif ()

if (PELEPHYSICS_SYCL)
    enable_language(SYCL)
    set(PELEPHYSICS_GPU_BACKEND SYCL)
endif ()

set(AMREX_HOME "" CACHE PATH "Path to AMReX top-level source directory")

if (AMREX_HOME)  # SUPERBUILD MODE

  message(STATUS "AMReX top-level source directory path: ${AMREX_HOME}")

  # Add AMReX module path
  list(APPEND CMAKE_MODULE_PATH ${AMREX_HOME}/Tools/CMake)
  set(AMREX_BINARY_DIR "${CMAKE_BINARY_DIR}/amrex" CACHE INTERNAL "Path to AMReX top-level binary directory")

  # Set required settings for AMReX
  set(USE_XSDK_DEFAULTS       ON                         CACHE INTERNAL "")
  set(XSDK_ENABLE_Fortran     OFF                        CACHE INTERNAL "")
  set(XSDK_PRECISION          ${PELEPHYSICS_PRECISION}   CACHE INTERNAL "")
  set(AMReX_SPACEDIM          ${PELEPHYSICS_SPACEDIM}    CACHE INTERNAL "")
  set(AMReX_MPI               ${PELEPHYSICS_MPI}         CACHE INTERNAL "")
  set(AMReX_OMP               ${PELEPHYSICS_OMP}         CACHE INTERNAL "")
  set(AMReX_GPU_BACKEND       ${PELEPHYSICS_GPU_BACKEND} CACHE INTERNAL "")
  set(AMReX_SUNDIALS          ON                         CACHE INTERNAL "")
  set(AMReX_EB                OFF                        CACHE INTERNAL "")
  set(AMReX_FORTRAN           OFF                        CACHE INTERNAL "")
  set(AMReX_LINEAR_SOLVERS    OFF                        CACHE INTERNAL "")
  set(AMReX_AMRLEVEL          OFF                        CACHE INTERNAL "")
  set(AMReX_BUILD_TUTORIALS   OFF                        CACHE INTERNAL "")

  # Add AMReX as a sub-project: incflo inherits all AMReX options
  add_subdirectory(${AMREX_HOME} ${AMREX_BINARY_DIR})

  message(STATUS "AMReX binary directory: ${AMREX_BINARY_DIR}")

else ()

#
# if AMReX::amrex is not part of the targets, use find_package
#
  if (NOT TARGET AMReX::amrex)
     set(AMREX_REQUIRED_COMPONENTS ${PELEPHYSICS_SPACEDIM}D ${PELEPHYSICS_PRECISION})
     if (PELEPHYSICS_MPI)
        list(APPEND AMREX_REQUIRED_COMPONENTS MPI)
     endif ()
     if (PELEPHYSICS_OMP)
        list(APPEND AMREX_REQUIRED_COMPONENTS OMP)
     endif ()
     if (NOT PELEPHYSICS_GPU_BACKEND STREQUAL "NONE")
         list(APPEND AMREX_REQUIRED_COMPONENTS ${PELEPHYSICS_GPU_BACKEND})
     endif ()
     find_package(AMReX CONFIG REQUIRED ${AMREX_REQUIRED_COMPONENTS} )
  endif ()

endif ()

set(SUNDIALS_HOME "" CACHE PATH "Path to SUNDIALS top-level source directory")

if (SUNDIALS_HOME)  # SUPERBUILD MODE
  message(STATUS "SUNDIALS top-level source directory path: ${SUNDIALS_HOME}")

  # Add SUNDIALS module path
  list(APPEND CMAKE_MODULE_PATH ${SUNDIALS_HOME}/cmake)
  set(SUNDIALS_BINARY_DIR "${CMAKE_BINARY_DIR}/sundials" CACHE INTERNAL "Path to SUNDIALS top-level binary directory")

  # Set required settings for SUNDIALS
  set(SUNDIALS_PRECISION   ${PELEPHYSICS_PRECISION} CACHE STRING "Floating point precision" FORCE)
  set(BUILD_ARKODE         ON                       CACHE INTERNAL "")
  set(BUILD_CVODE          ON                       CACHE INTERNAL "")
  set(SUNDIALS_INDEX_SIZE  32                       CACHE INTERNAL "")
  set(BUILD_CVODES         OFF                      CACHE INTERNAL "")
  set(BUILD_IDA            OFF                      CACHE INTERNAL "")
  set(BUILD_IDAS           OFF                      CACHE INTERNAL "")
  set(BUILD_KINSOL         OFF                      CACHE INTERNAL "")
  set(BUILD_CPODES         OFF                      CACHE INTERNAL "")
  set(BUILD_EXAMPLES       OFF                      CACHE INTERNAL "")
  set(_BUILD_EXAMPLES      OFF                      CACHE INTERNAL "")
  set(BUILD_TESTING        OFF                      CACHE INTERNAL "")
  set(EXAMPLES_ENABLE_C    OFF                      CACHE INTERNAL "")
  set(EXAMPLES_ENABLE_CXX  OFF                      CACHE INTERNAL "")
  set(EXAMPLES_INSTALL     OFF                      CACHE INTERNAL "")
  set(ENABLE_CUDA          ${PELEPHYSICS_CUDA}      CACHE INTERNAL "")
  set(ENABLE_HIP           ${PELEPHYSICS_HIP}       CACHE INTERNAL "")
  set(ENABLE_SYCL          ${PELEPHYSICS_SYCL}      CACHE INTERNAL "")
  set(CMAKE_POLICY_DEFAULT_CMP0077 NEW)

  add_subdirectory(${SUNDIALS_HOME} ${SUNDIALS_BINARY_DIR})

  message(STATUS "SUNDIALS binary directory: ${SUNDIALS_BINARY_DIR}")
else ()

#
# Find SUNDIALS, check for required components
#
  if (NOT TARGET SUNDIALS::cvode)
    # Check if we can require components
    find_package(SUNDIALS)
  endif ()

endif ()

#
# Enable CUDA if requested
#
if (PELEPHYSICS_GPU_BACKEND STREQUAL "CUDA")
   include(AMReXTargetHelpers)
endif ()

#
# Define the object library to compile
#
set(pelephysics_lib_name PelePhysics-${PELEPHYSICS_CHEMISTRY_MODEL}-${PELEPHYSICS_EOS_MODEL}-${PELEPHYSICS_TRANSPORT_MODEL})
add_library(${pelephysics_lib_name} OBJECT)
target_link_libraries(${pelephysics_lib_name} PUBLIC AMReX::amrex)
target_link_libraries(${pelephysics_lib_name} PUBLIC SUNDIALS::cvode SUNDIALS::arkode)

add_subdirectory(Source)

set(PELEPHYSICS_MECH_DIR "Support/Mechanism/")
add_subdirectory(${PELEPHYSICS_MECH_DIR})

if (PELEPHYSICS_UTILS)
  add_subdirectory(Utility/PMF)
  add_subdirectory(Utility/TurbInflow)
  add_subdirectory(Utility/PltFileManager)
  add_subdirectory(Utility/Diagnostics)
endif ()

if (PELEPHYSICS_TRANSPORT)
  add_subdirectory(Transport)
endif ()

if (PELEPHYSICS_EOS)
  add_subdirectory(Eos)
endif ()

if (PELEPHYSICS_REACTIONS)
  add_subdirectory(Reactions)
endif ()

if (PELEPHYSICS_GPU_BACKEND STREQUAL "CUDA")
  setup_target_for_cuda_compilation(${pelephysics_lib_name})
endif ()

if ( NOT CMAKE_CXX_FLAGS )
  target_link_libraries(${pelephysics_lib_name} PUBLIC AMReX::Flags_CXX)
endif ()

# Installation rules
include(CMakePackageConfigHelpers)
include(GNUInstallDirs)

message( STATUS "PelePhysics configuration summary: ")
message( STATUS "   Build type               = ${CMAKE_BUILD_TYPE}")
message( STATUS "   Install directory        = ${CMAKE_INSTALL_PREFIX}")
message( STATUS "   CMake target             = ${pelephysics_lib_name}")

# Collect all headers and make them installable with the target
get_target_property(PELEPHYSICS_INCLUDES ${pelephysics_lib_name} SOURCES)
list(FILTER PELEPHYSICS_INCLUDES INCLUDE REGEX "\\.H")
set_target_properties(
    ${pelephysics_lib_name} PROPERTIES PUBLIC_HEADER "${PELEPHYSICS_INCLUDES}")

# Install PelePhysics
install(
  TARGETS ${pelephysics_lib_name}
  EXPORT ${PROJECT_NAME}Targets
  RUNTIME       DESTINATION bin
  ARCHIVE       DESTINATION lib
  LIBRARY       DESTINATION lib
  INCLUDES      DESTINATION include
  PUBLIC_HEADER DESTINATION include
  )

# Make PelePhysics discoverable using `find_package`
install(
  EXPORT ${PROJECT_NAME}Targets
  NAMESPACE ${PROJECT_NAME}::
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  )
configure_package_config_file(
  CMake/${PROJECT_NAME}Config.cmake.in
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  INSTALL_DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  )
install(FILES
  ${CMAKE_CURRENT_BINARY_DIR}/${PROJECT_NAME}Config.cmake
  DESTINATION ${CMAKE_INSTALL_LIBDIR}/cmake/${PROJECT_NAME}
  )
