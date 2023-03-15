target_include_directories(${pelephysics_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(${pelephysics_lib_name}
  PRIVATE
    EOS.H
    EOS.cpp
    Fuego.H
    GammaLaw.H
    SRK.H
  )

if("${PELEPHYSICS_EOS_MODEL}" STREQUAL "GammaLaw")
  target_compile_definitions(${pelephysics_lib_name} PUBLIC USE_GAMMALAW_EOS)
endif()
if("${PELEPHYSICS_EOS_MODEL}" STREQUAL "Fuego")
  target_compile_definitions(${pelephysics_lib_name} PUBLIC USE_FUEGO_EOS)
endif()
if("${PELEPHYSICS_EOS_MODEL}" STREQUAL "Soave-Redlich-Kwong")
  target_compile_definitions(${pelephysics_lib_name} PUBLIC USE_SRK_EOS)
endif()
