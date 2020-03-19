target_sources(${pele_exe_name}
   PRIVATE
      mechanism.h
      mechanism.cpp
   )

target_include_directories(${pele_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
