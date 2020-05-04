# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_LOKI
    DESCRIPTION "Use Loki/CLAW source-to-source transformations with CLOUDSC " DEFAULT OFF
    CONDITION Serialbox_FOUND
)

if( HAVE_CLOUDSC_LOKI )

    # Define the binary build target for this variant
    ecbuild_add_executable( TARGET dwarf-cloudsc-loki
        SOURCES
            dwarf_cloudsc.F90
            cloudsc_global_state_mod.F90
            cloudsc_driver_loki_mod.F90
            cloudsc.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries( dwarf-cloudsc-loki PRIVATE cloudsc-common-lib )

    if( TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc-loki PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )
endif()
