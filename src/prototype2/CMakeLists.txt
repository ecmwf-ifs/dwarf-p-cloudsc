ecbuild_enable_fortran(REQUIRED MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

find_package( OpenMP COMPONENTS Fortran )
find_package( Serialbox )

if( Serialbox_FOUND AND SERIALBOX_HAS_FORTRAN )

    ecbuild_add_library(
        TARGET  dwarf-cloudsc-v2-lib
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
                           Serialbox::Serialbox_Fortran
                           OpenMP::OpenMP_Fortran )

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-v2
        SOURCES
            dwarf_cloudsc.F90 
    )
    target_link_libraries( dwarf-cloudsc-v2 PRIVATE
                           dwarf-cloudsc-v2-lib
                           OpenMP::OpenMP_Fortran )


else()
    message(STATUS "Serialbox not found, disabling prototype2" )
endif()
