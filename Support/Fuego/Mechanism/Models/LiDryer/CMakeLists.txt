target_sources(${pele_exe_name}
   PRIVATE
      chemistry_file.H
      mechanism.h
      LiDryer.cpp
   )

target_include_directories(${pele_exe_name} PRIVATE ${CMAKE_CURRENT_SOURCE_DIR})
