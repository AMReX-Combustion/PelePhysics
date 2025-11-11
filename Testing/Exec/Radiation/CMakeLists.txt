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

if(PELE_ENABLE_MPI)
  set(MPI_COMMANDS "${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} 2 ${MPIEXEC_PREFLAGS}")
endif()

if(PELE_ENABLE_EB)
  add_executable(PeleRad_POneSingleEB.exe tstPOneSingleEB.cpp)
  target_link_libraries(PeleRad_POneSingleEB.exe PRIVATE PeleRad)
  add_test(NAME PeleRad_POneSingleEB_Test COMMAND bash -c "${MPI_COMMANDS} ${CMAKE_CURRENT_BINARY_DIR}/PeleRad_POneSingleEB.exe ${MPIEXEC_POSTFLAGS} ${CMAKE_CURRENT_SOURCE_DIR}/inputs/inputs.tstPOneSingleEB")
  set_tests_properties(PeleRad_POneSingleEB_Test PROPERTIES TIMEOUT 7200)
  
  add_executable(PeleRad_POneMultiEB.exe tstPOneMultiEB.cpp)
  target_link_libraries(PeleRad_POneMultiEB.exe PRIVATE PeleRad)
  add_test(NAME PeleRad_POneMultiEB_Test COMMAND bash -c "${MPI_COMMANDS} ${CMAKE_CURRENT_BINARY_DIR}/PeleRad_POneMultiEB.exe ${MPIEXEC_POSTFLAGS} ${CMAKE_CURRENT_SOURCE_DIR}/inputs/inputs.tstPOneMultiEB")
  set_tests_properties(PeleRad_POneMultiEB_Test PROPERTIES TIMEOUT 7200)
else()
  add_executable(PeleRad_POneSingle.exe tstPOneSingle.cpp)
  target_link_libraries(PeleRad_POneSingle.exe PRIVATE PeleRad AMReX::amrex)
  add_test(NAME PeleRad_POneSingle_Test COMMAND bash -c "${MPI_COMMANDS} ${CMAKE_CURRENT_BINARY_DIR}/PeleRad_POneSingle.exe ${MPIEXEC_POSTFLAGS} ${CMAKE_CURRENT_SOURCE_DIR}/inputs/inputs.tstPOneSingle")
  set_tests_properties(PeleRad_POneSingle_Test PROPERTIES TIMEOUT 7200)
 
  add_executable(PeleRad_POneSingleAF.exe tstPOneSingleAF.cpp)
  target_link_libraries(PeleRad_POneSingleAF.exe PRIVATE PeleRad)
  add_test(NAME PeleRad_POneSingleAF_Test COMMAND bash -c "${MPI_COMMANDS} ${CMAKE_CURRENT_BINARY_DIR}/PeleRad_POneSingleAF.exe ${MPIEXEC_POSTFLAGS} ${CMAKE_CURRENT_SOURCE_DIR}/inputs/inputs.tstPOneSingleAF")
  set_tests_properties(PeleRad_POneSingleAF_Test PROPERTIES TIMEOUT 7200)
  
  add_executable(PeleRad_POneMultiAF.exe tstPOneMultiAF.cpp)
  target_link_libraries(PeleRad_POneMultiAF.exe PRIVATE PeleRad)
  add_test(NAME PeleRad_POneMultiAF_Test_Composite COMMAND bash -c "${MPI_COMMANDS} ${CMAKE_CURRENT_BINARY_DIR}/PeleRad_POneMultiAF.exe ${MPIEXEC_POSTFLAGS} ${CMAKE_CURRENT_SOURCE_DIR}/inputs/inputs.tstPOneMultiAF")
  set_tests_properties(PeleRad_POneMultiAF_Test_Composite PROPERTIES TIMEOUT 7200)
  add_test(NAME PeleRad_POneMultiAF_Test_LevbyLev COMMAND bash -c "${MPI_COMMANDS} ${CMAKE_CURRENT_BINARY_DIR}/PeleRad_POneMultiAF.exe ${MPIEXEC_POSTFLAGS} ${CMAKE_CURRENT_SOURCE_DIR}/inputs/inputs.tstPOneMultiAFLevbyLev")
  set_tests_properties(PeleRad_POneMultiAF_Test_LevbyLev PROPERTIES TIMEOUT 7200)

  if(PELE_DIM GREATER 2)
    add_executable(PeleRad_POneMulti.exe tstPOneMulti.cpp)
    target_link_libraries(PeleRad_POneMulti.exe PRIVATE PeleRad)
    add_test(NAME PeleRad_POneMulti_Test_Composite COMMAND bash -c "${MPI_COMMANDS} ${CMAKE_CURRENT_BINARY_DIR}/PeleRad_POneMulti.exe ${MPIEXEC_POSTFLAGS} ${CMAKE_CURRENT_SOURCE_DIR}/inputs/inputs.tstPOneMulti")
    set_tests_properties(PeleRad_POneMulti_Test_Composite PROPERTIES TIMEOUT 7200)
    add_test(NAME PeleRad_POneMulti_Test_LevbyLev COMMAND bash -c "${MPI_COMMANDS} ${CMAKE_CURRENT_BINARY_DIR}/PeleRad_POneMulti.exe ${MPIEXEC_POSTFLAGS} ${CMAKE_CURRENT_SOURCE_DIR}/inputs/inputs.tstPOneMultiLevbyLev")
    set_tests_properties(PeleRad_POneMulti_Test_LevbyLev PROPERTIES TIMEOUT 7200)
  endif()
endif()

if(PELE_ENABLE_CUDA)
  if(PELE_ENABLE_EB)
    setup_target_for_cuda_compilation(PeleRad_POneSingleEB.exe)
    setup_target_for_cuda_compilation(PeleRad_POneMultiEB.exe)
  else()
    setup_target_for_cuda_compilation(PeleRad_POneSingle.exe)
    setup_target_for_cuda_compilation(PeleRad_POneMulti.exe)
    setup_target_for_cuda_compilation(PeleRad_POneSingleAF.exe)
    setup_target_for_cuda_compilation(PeleRad_POneMultiAF.exe)
  endif()
endif()
