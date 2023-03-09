target_include_directories(${pelephysics_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(${pelephysics_lib_name}
  PRIVATE
    PMF.H
    PMFData.cpp
    PMFData.H
  )
