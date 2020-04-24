ecbuild_enable_fortran(REQUIRED MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module)
include_directories(${CMAKE_CURRENT_SOURCE_DIR}/include)

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_FORTRAN
    DESCRIPTION "Build the new Fortran version CLOUDSC using Serialbox" DEFAULT ON
)

if( HAVE_CLOUDSC_FORTRAN )

    # Define the binary build target for this variant
    ecbuild_add_executable( TARGET dwarf-cloudsc-fortran
        SOURCES
            dwarf_cloudsc.F90
            cloudsc_global_state_mod.F90
            cloudsc_driver_mod.F90
            cloudsc.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries( dwarf-cloudsc-fortran PRIVATE cloudsc-common-lib )

    if( OpenMP_Fortran_FOUND )
        target_link_libraries( dwarf-cloudsc-fortran PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )
endif()
