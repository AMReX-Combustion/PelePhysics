cmake_minimum_required (VERSION 3.10)
project(PelePhysics-Docs NONE)
list(APPEND CMAKE_MODULE_PATH "${CMAKE_SOURCE_DIR}/cmake")

find_package(Sphinx REQUIRED)

add_custom_target(sphinx ALL
  COMMAND ${SPHINX_EXECUTABLE} -M html
  "${CMAKE_CURRENT_SOURCE_DIR}/sphinx" "${CMAKE_CURRENT_BINARY_DIR}/sphinx"
  -c ${CMAKE_CURRENT_SOURCE_DIR}/sphinx -W --keep-going -n)
