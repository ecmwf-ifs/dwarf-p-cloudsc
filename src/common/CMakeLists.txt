ecbuild_enable_fortran(REQUIRED MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package( OpenMP COMPONENTS Fortran )
find_package( Serialbox )

# Add option for single-precision builds
ecbuild_add_option(
    FEATURE SINGLE_PRECISION
    DESCRIPTION "Build CLOUDSC in single precision" DEFAULT OFF
)
if( HAVE_SINGLE_PRECISION )
  list(APPEND CLOUDSC_DEFINITIONS SINGLE)
endif()


if( Serialbox_FOUND AND SERIALBOX_HAS_FORTRAN )

    ecbuild_add_library(
        TARGET cloudsc-common-lib
        SOURCES 
            module/parkind1.F90
            module/yoecldp.F90
            module/yomcst.F90
            module/yoethf.F90
            module/yoephli.F90
            module/yomphyder.F90
            module/routines.F90
            module/abor1.F90
            module/timer_mod.F90
            module/mycpu.c
            module/ec_pmon_mod.F90
            module/load_array_mod.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries(
        cloudsc-common-lib
        PRIVATE
            Serialbox::Serialbox_Fortran
            OpenMP::OpenMP_Fortran
    )

else()
    message(STATUS "Serialbox not found, disabling cloudsc-common-lib" )
endif()
