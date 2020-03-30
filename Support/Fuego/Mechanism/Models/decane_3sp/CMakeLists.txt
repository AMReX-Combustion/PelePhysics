target_sources(${pele_exe_name}
   PRIVATE
      chemistry_file.h
      mechanism.h
      decane_3sp.cpp
   )

target_include_directories(${pele_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
