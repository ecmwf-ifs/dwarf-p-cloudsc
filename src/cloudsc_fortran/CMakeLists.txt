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

    ecbuild_add_executable(
        TARGET dwarf-cloudsc-fortran
        SOURCES
            dwarf_cloudsc.F90
            cloudsc_global_state_mod.F90
            cloudsc_driver_mod.F90
            cloudsc.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries(
        dwarf-cloudsc-fortran
	PRIVATE
            cloudsc-common-lib
            OpenMP::OpenMP_Fortran
    )

    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )

else()
    message(STATUS "Serialbox not found, disabling prototype2" )
endif()
