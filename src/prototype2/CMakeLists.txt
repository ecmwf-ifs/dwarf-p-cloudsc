ecbuild_enable_fortran(REQUIRED MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Note, that we link statically here, since Serialbox only exposes static Fortran libs.
# Serailbox is used to generate the off-line input and reference data for alternative
# versions of this dwarf, eg. prototype2.
#
# We should return to dynamic linking once Serialbox::SerialboxFortranShared is available.

if( Serialbox_FOUND )
ecbuild_add_library(
    TARGET  dwarf-cloudsc-v2-lib
    TYPE STATIC
    DEFINITIONS
        ${CLOUDSC_DEFINITIONS}
    SOURCES 
        modules/parkind1.F90
        modules/yoecldp.F90
        modules/yomcst.F90
        modules/yoethf.F90
        modules/yoephli.F90
        modules/yomphyder.F90
        modules/routines.F90
        modules/abor1.F90
        modules/timer_mod.F90
        modules/mycpu.c
        cloudsc/load_state_mod.F90
        cloudsc/cloudsc_global_state_mod.F90
        cloudsc/cloudsc_driver_mod.F90
        cloudsc/cloudsc.F90
)
target_link_libraries( dwarf-cloudsc-v2-lib PRIVATE
                       Serialbox::SerialboxFortranStatic
                       OpenMP::OpenMP_Fortran )

ecbuild_add_executable(
    TARGET  dwarf-cloudsc-v2
    SOURCES
        dwarf_cloudsc.F90 
    LIBS
        dwarf-cloudsc-v2-lib
)
else()
  message(STATUS "Serialbox not found, disabling prototype2" )
endif()
