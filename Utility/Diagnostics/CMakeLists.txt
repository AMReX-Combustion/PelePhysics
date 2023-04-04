target_include_directories(${pelephysics_lib_name} PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_LIST_DIR}>)

target_sources(${pelephysics_lib_name}
    PRIVATE
        DiagBase.cpp
        DiagBase.H
        DiagConditional.cpp
        DiagConditional.H
        DiagFilter.cpp
        DiagFilter.H
        DiagFramePlane.cpp
        DiagFramePlane.H
        DiagPDF.cpp
        DiagPDF.H
    )
