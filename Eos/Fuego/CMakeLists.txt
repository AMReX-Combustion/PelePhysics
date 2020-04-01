target_sources(${pele_exe_name}
   PRIVATE
      EOS.cpp
      EOS.H
   )

target_include_directories(${pele_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
