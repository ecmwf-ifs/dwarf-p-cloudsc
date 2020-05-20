# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_LOKI
    DESCRIPTION "Use Loki/CLAW source-to-source transformations with CLOUDSC " DEFAULT OFF
    CONDITION Serialbox_FOUND
)

if( HAVE_CLOUDSC_LOKI )

    ####################################################
    ##  Define various pre-processing modes via Loki  ##
    ####################################################

    set( COMMON_MODULE "${CMAKE_CURRENT_SOURCE_DIR}/../common/module" )
    set( COMMON_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/../common/include" )
    set( XMOD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/xmod" )

    set( LOKI_FRONTEND "ofp" CACHE STRING "Frontend parser for Loki transforms" )

    # Ensure xmod directory for OMNI frontend
    file(MAKE_DIRECTORY ${XMOD_DIR})

    # "Idempotence mode" for debugging Loki internals
    add_custom_command(
        OUTPUT cloudsc.idem.F90 cloudsc_driver_loki_mod.idem.F90
        COMMAND loki-transform.py idem --frontend ${LOKI_FRONTEND}
                      --header ${COMMON_MODULE}/yomphyder.F90
                      --driver ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_driver_loki_mod.F90
                      --source ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc.F90
                      --include ${COMMON_INCLUDE}
                      --xmod ${XMOD_DIR}
                      --out-path ${CMAKE_CURRENT_BINARY_DIR}
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90
	COMMENT "[Loki] Pre-processing: mode=idem frontend=${LOKI_FRONTEND}"
    )

    ####################################################
    ##  Define build targets for Loki transformations ##
    ####################################################

    # Define the binary build target for this variant
    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-idem
        SOURCES
            dwarf_cloudsc.F90
            cloudsc_global_state_mod.F90
            cloudsc_driver_loki_mod.idem.F90
            cloudsc.idem.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries( dwarf-cloudsc-loki-idem PRIVATE cloudsc-common-lib )

    if( TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc-loki-idem PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )
endif()
