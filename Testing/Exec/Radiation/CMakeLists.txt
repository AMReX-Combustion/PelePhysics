add_library(PeleRad OBJECT)

set(RADIATION_INCLUDE ${CMAKE_SOURCE_DIR}/Submodules/PelePhysics/Source/Radiation)
target_sources(PeleRad
  PRIVATE
    ${RADIATION_INCLUDE}/AMRParam.H
    ${RADIATION_INCLUDE}/Constants.H
    ${RADIATION_INCLUDE}/MLMGParam.H
    ${RADIATION_INCLUDE}/POneMulti.H
    ${RADIATION_INCLUDE}/POneMultiEB.H
    ${RADIATION_INCLUDE}/POneMultiLevbyLev.H
    ${RADIATION_INCLUDE}/POneSingle.H
    ${RADIATION_INCLUDE}/POneSingleEB.H
    ${RADIATION_INCLUDE}/PeleCRad.H
    ${RADIATION_INCLUDE}/PeleLMRad.H
    ${RADIATION_INCLUDE}/PlanckMean.H
    ${RADIATION_INCLUDE}/SpectralModels.H
)
target_include_directories(PeleRad PUBLIC ${RADIATION_INCLUDE})
target_link_libraries(PeleRad PUBLIC AMReX::amrex)

file(COPY kpDB DESTINATION "${CMAKE_CURRENT_BINARY_DIR}/")

if(PELELMEX_ENABLE_EB)
  add_executable(PeleRad_POneSingleEB.exe tstPOneSingleEB.cpp)
  target_link_libraries(PeleRad_POneSingleEB.exe PRIVATE PeleRad)
  add_test(NAME PeleRad_POneSingleEB_Test COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/PeleRad_POneSingleEB.exe ${CMAKE_CURRENT_SOURCE_DIR}/inputs/inputs.tstPOneSingleEB)
  
  add_executable(PeleRad_POneMultiEB.exe tstPOneMultiEB.cpp)
  target_link_libraries(PeleRad_POneMultiEB.exe PRIVATE PeleRad)
  add_test(NAME PeleRad_POneMultiEB_Test COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/PeleRad_POneMultiEB.exe ${CMAKE_CURRENT_SOURCE_DIR}/inputs/inputs.tstPOneMultiEB)
else()
  add_executable(PeleRad_POneSingle.exe tstPOneSingle.cpp)
  target_link_libraries(PeleRad_POneSingle.exe PRIVATE PeleRad AMReX::amrex)
  add_test(NAME PeleRad_POneSingle_Test COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/PeleRad_POneSingle.exe ${CMAKE_CURRENT_SOURCE_DIR}/inputs/inputs.tstPOneSingle)
    
  add_executable(PeleRad_POneMulti.exe tstPOneMulti.cpp)
  target_link_libraries(PeleRad_POneMulti.exe PRIVATE PeleRad)
  add_test(NAME PeleRad_POneMulti_Test_Composite COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/PeleRad_POneMulti.exe ${CMAKE_CURRENT_SOURCE_DIR}/inputs/inputs.tstPOneMulti)
  add_test(NAME PeleRad_POneMulti_Test_LevbyLev COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/PeleRad_POneMulti.exe ${CMAKE_CURRENT_SOURCE_DIR}/inputs/inputs.tstPOneMultiLevbyLev)
  
  add_executable(PeleRad_POneSingleAF.exe tstPOneSingleAF.cpp)
  target_link_libraries(PeleRad_POneSingleAF.exe PRIVATE PeleRad)
  add_test(NAME PeleRad_POneSingleAF_Test COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/PeleRad_POneSingleAF.exe ${CMAKE_CURRENT_SOURCE_DIR}/inputs/inputs.tstPOneSingleAF)
  
  add_executable(PeleRad_POneMultiAF.exe tstPOneMultiAF.cpp)
  target_link_libraries(PeleRad_POneMultiAF.exe PRIVATE PeleRad)
  add_test(NAME PeleRad_POneMultiAF_Test_Composite COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/PeleRad_POneMultiAF.exe ${CMAKE_CURRENT_SOURCE_DIR}/inputs/inputs.tstPOneMultiAF)
  add_test(NAME PeleRad_POneMultiAF_Test_LevbyLev COMMAND ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} ${MPIEXEC_MAX_NUMPROCS} ${MPIEXEC_PREFLAGS} ${CMAKE_CURRENT_BINARY_DIR}/PeleRad_POneMultiAF.exe ${CMAKE_CURRENT_SOURCE_DIR}/inputs/inputs.tstPOneMultiAFLevbyLev)
endif()

if(PELELMEX_ENABLE_CUDA)
  setup_target_for_cuda_compilation(PeleRad_POneSingle.exe)
  setup_target_for_cuda_compilation(PeleRad_POneSingleEB.exe)
  setup_target_for_cuda_compilation(PeleRad_POneMultiEB.exe)
  setup_target_for_cuda_compilation(PeleRad_POneMulti.exe)
  setup_target_for_cuda_compilation(PeleRad_POneSingleAF.exe)
  setup_target_for_cuda_compilation(PeleRad_POneMultiAF.exe)
endif()
