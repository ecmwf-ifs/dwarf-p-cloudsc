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
        OUTPUT loki-idem/cloudsc.idem.F90 loki-idem/cloudsc_driver_loki_mod.idem.F90
        COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/loki-idem
        COMMAND loki-transform.py idem --frontend ${LOKI_FRONTEND}
                      --header ${COMMON_MODULE}/yomphyder.F90
                      --driver ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_driver_loki_mod.F90
                      --source ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc.F90
                      --include ${COMMON_INCLUDE}
                      --xmod ${XMOD_DIR}
                      --out-path ${CMAKE_CURRENT_BINARY_DIR}/loki-idem
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90
	COMMENT "[Loki] Pre-processing: mode=idem frontend=${LOKI_FRONTEND}"
    )

    # "Single Column Abstraction (SCA) mode" for removing horizontal vector loops
    add_custom_command(
        OUTPUT loki-sca/cloudsc.sca.F90 loki-sca/cloudsc_driver_loki_mod.sca.F90
        COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/loki-sca
        COMMAND loki-transform.py convert --mode sca --frontend ${LOKI_FRONTEND}
                      --header ${COMMON_MODULE}/yomphyder.F90
                      --driver ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_driver_loki_mod.F90
                      --source ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc.F90
                      --include ${COMMON_INCLUDE}
                      --xmod ${XMOD_DIR}
                      --out-path ${CMAKE_CURRENT_BINARY_DIR}/loki-sca
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90
	COMMENT "[Loki] Pre-processing: mode=sca frontend=${LOKI_FRONTEND}"
    )

    # CLAW-GPU mode for generating OpenACC-parallel GPU code via loki-generated SCA
    add_custom_command(
        OUTPUT loki-claw-gpu/cloudsc.claw.gpu.F90 loki-claw-gpu/cloudsc_driver_loki_mod.claw.gpu.F90
        COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/loki-claw-gpu
        COMMAND loki-transform.py convert --mode claw --frontend ${LOKI_FRONTEND}
                      --header ${COMMON_MODULE}/yomphyder.F90
                      --driver ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_driver_loki_mod.F90
                      --source ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc.F90
                      --include ${COMMON_INCLUDE}
                      --xmod ${XMOD_DIR}
                      --out-path ${CMAKE_CURRENT_BINARY_DIR}/loki-claw-gpu
        COMMAND clawfc -w 132 --target=gpu --directive=openacc
	              --model-config=${CMAKE_CURRENT_SOURCE_DIR}/claw_cloudsc.config
                      -I ${COMMON_INCLUDE} -J ${XMOD_DIR}
                      -o loki-claw-gpu/cloudsc.claw.gpu.F90
                      loki-claw-gpu/cloudsc.claw.F90
        COMMAND clawfc -w 132 --target=gpu --directive=openacc
	              --model-config=${CMAKE_CURRENT_SOURCE_DIR}/claw_cloudsc.config
                      -I ${COMMON_INCLUDE} -J ${XMOD_DIR}
                      -o loki-claw-gpu/cloudsc_driver_loki_mod.claw.gpu.F90
                      loki-claw-gpu/cloudsc_driver_loki_mod.claw.F90
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90
	COMMENT "[Loki] Pre-processing: mode=claw-gpu frontend=${LOKI_FRONTEND}"
    )


    # CLAW-CPU mode for generating vectorized CPU host code via loki-generated SCA
    add_custom_command(
        OUTPUT loki-claw-cpu/cloudsc.claw.cpu.F90 loki-claw-cpu/cloudsc_driver_loki_mod.claw.cpu.F90
        COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/loki-claw-cpu
        COMMAND loki-transform.py convert --mode claw --frontend ${LOKI_FRONTEND}
                      --header ${COMMON_MODULE}/yomphyder.F90
                      --driver ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_driver_loki_mod.F90
                      --source ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc.F90
                      --include ${COMMON_INCLUDE}
                      --xmod ${XMOD_DIR}
                      --out-path ${CMAKE_CURRENT_BINARY_DIR}/loki-claw-cpu
        # We purposefully suppress CLAWs insertion of OpenMP loops,
	# as their are already present in the outer driver.
        COMMAND clawfc -w 132 --target=cpu --directive=none
	              --model-config=${CMAKE_CURRENT_SOURCE_DIR}/claw_cloudsc.config
                      -I ${COMMON_INCLUDE} -J ${XMOD_DIR}
                      -o loki-claw-cpu/cloudsc.claw.cpu.F90
                      loki-claw-cpu/cloudsc.claw.F90
        COMMAND clawfc -w 132 --target=cpu --directive=none
	              --model-config=${CMAKE_CURRENT_SOURCE_DIR}/claw_cloudsc.config
                      -I ${COMMON_INCLUDE} -J ${XMOD_DIR}
                      -o loki-claw-cpu/cloudsc_driver_loki_mod.claw.cpu.F90
                      loki-claw-cpu/cloudsc_driver_loki_mod.claw.F90
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90
	COMMENT "[Loki] Pre-processing: mode=claw-cpu frontend=${LOKI_FRONTEND}"
    )


    ####################################################
    ##  Define build targets for Loki transformations ##
    ####################################################

    # Define the "idem" build target for this variant
    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-idem
        SOURCES
            dwarf_cloudsc.F90
            cloudsc_global_state_mod.F90
            loki-idem/cloudsc_driver_loki_mod.idem.F90
            loki-idem/cloudsc.idem.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries( dwarf-cloudsc-loki-idem PRIVATE cloudsc-common-lib )
    if( TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc-loki-idem PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    # Define the "sca" build target for this variant
    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-sca
        SOURCES
            dwarf_cloudsc.F90
            cloudsc_global_state_mod.F90
            loki-sca/cloudsc_driver_loki_mod.sca.F90
            loki-sca/cloudsc.sca.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries( dwarf-cloudsc-loki-sca PRIVATE cloudsc-common-lib )
    if( TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc-loki-sca PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    # Define the CLAW-GPU build target for this variant
    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-claw-gpu
        SOURCES
            dwarf_cloudsc.F90
            cloudsc_global_state_mod.F90
            loki-claw-gpu/cloudsc_driver_loki_mod.claw.gpu.F90
            loki-claw-gpu/cloudsc.claw.gpu.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries( dwarf-cloudsc-loki-claw-gpu PRIVATE cloudsc-common-lib )
    if( TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc-loki-claw-gpu PRIVATE OpenMP::OpenMP_Fortran )
    endif()
    if( TARGET OpenACC::OpenACC_Fortran )
        target_link_libraries( dwarf-cloudsc-loki-claw-gpu PRIVATE OpenACC::OpenACC_Fortran )
    endif()

    # Define the CLAW-CPU build target for this variant
    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-claw-cpu
        SOURCES
            dwarf_cloudsc.F90
            cloudsc_global_state_mod.F90
            loki-claw-cpu/cloudsc_driver_loki_mod.claw.cpu.F90
            loki-claw-cpu/cloudsc.claw.cpu.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries( dwarf-cloudsc-loki-claw-cpu PRIVATE cloudsc-common-lib )
    if( TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc-loki-claw-cpu PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    ##############################################################################
    # C-transpilation mode for generating vectorized C host code (experimental!)
    ##############################################################################
    add_custom_command(
	OUTPUT loki-c/cloudsc_driver_loki_mod.c.F90
               loki-c/cloudsc_fc.f90 loki-c/cloudsc_c.c
               loki-c/yoethf_fc.f90 loki-c/yomcst_fc.f90
               loki-c/yoecldp_fc.f90
        COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/loki-c
        COMMAND loki-transform.py transpile
                      --header ${COMMON_MODULE}/yomphyder.F90
                      --header ${COMMON_MODULE}/yomcst.F90
                      --header ${COMMON_MODULE}/yoethf.F90
                      --header ${COMMON_MODULE}/yoecldp.F90
                      --driver ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_driver_loki_mod.F90
                      --source ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc.F90
                      --include ${COMMON_INCLUDE}
                      --xmod ${XMOD_DIR}
                      --out-path ${CMAKE_CURRENT_BINARY_DIR}/loki-c
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90
	COMMENT "[Loki] Pre-processing: mode=transpile frontend=omni"
    )

    # Define the CLAW-CPU build target for this variant
    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-c
        SOURCES
            dwarf_cloudsc.F90
            cloudsc_global_state_mod.F90
            loki-c/cloudsc_driver_loki_mod.c.F90
            loki-c/cloudsc_fc.f90
            loki-c/cloudsc_c.c
            loki-c/yoethf_fc.f90
            loki-c/yomcst_fc.f90
            loki-c/yoecldp_fc.f90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries( dwarf-cloudsc-loki-c PRIVATE cloudsc-common-lib )
    if( TARGET OpenMP::OpenMP_Fortran AND TARGET OpenMP::OpenMP_C )
        target_link_libraries( dwarf-cloudsc-loki-c PRIVATE OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc-loki-c PRIVATE OpenMP::OpenMP_C )
    endif()

    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )
endif()
