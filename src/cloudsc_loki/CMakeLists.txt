# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_LOKI
    DESCRIPTION "Use Loki source-to-source transformations with CLOUDSC " DEFAULT OFF
    CONDITION Serialbox_FOUND OR HDF5_FOUND
)
ecbuild_add_option( FEATURE CLOUDSC_LOKI_CLAW
    DESCRIPTION "Use CLAW-based Loki source-to-source transformations" DEFAULT ON
    CONDITION HAVE_CLOUDSC_LOKI
)

if( HAVE_CLOUDSC_LOKI )

    ####################################################
    ##  Define various pre-processing modes via Loki  ##
    ####################################################

    set( COMMON_MODULE "${CMAKE_CURRENT_SOURCE_DIR}/../common/module" )
    set( COMMON_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/../common/include" )
    set( XMOD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/xmod" )

    set( LOKI_FRONTEND "fp" CACHE STRING "Frontend parser for Loki transforms" )

    # OFP frontend cannot deal with statement functions, so we toggle them here
    set( CLOUDSC_DEFINE_STMT_FUNC "" )
    if(NOT "${LOKI_FRONTEND}" STREQUAL "ofp")
      set( CLOUDSC_DEFINE_STMT_FUNC --define CLOUDSC_STMT_FUNC )
    endif()

    # Ensure xmod directory for OMNI frontend
    file(MAKE_DIRECTORY ${XMOD_DIR})

    if ( TARGET clawfc AND ${LOKI_FRONTEND} STREQUAL "omni" )
        # Ugly hack but I don't have a better solution: We need to add F_FRONT
        # (which is installed in the same directory as clawfc) to the PATH, if
        # OMNI is used as a frontend. Hence we have to update the environment in the below
        # add_custom_command calls to loki-transform.py.
        # Unfortunately, this environment updated breaks the CMake feature of recognizing
        # the executable in add_custom_command and choosing the correct path if it was
        # previously declared as an executable. Therefore, we have to insert also
        # loki-transform.py into the PATH variable.
        get_target_property( _LOKI_TRANSFORM_EXECUTABLE loki-transform.py IMPORTED_LOCATION )
        get_filename_component( _LOKI_TRANSFORM_LOCATION ${_LOKI_TRANSFORM_EXECUTABLE} DIRECTORY )
        get_target_property( _CLAWFC_EXECUTABLE clawfc IMPORTED_LOCATION )
        get_filename_component( _CLAWFC_LOCATION ${_CLAWFC_EXECUTABLE} DIRECTORY )
        set( _LOKI_TRANSFORM_ENV PATH=${_LOKI_TRANSFORM_LOCATION}:${_CLAWFC_LOCATION}:$ENV{PATH} )
        set( _LOKI_TRANSFORM ${CMAKE_COMMAND} -E env ${_LOKI_TRANSFORM_ENV} loki-transform.py )
    else()
        # This is how it is meant to be: We can rely on CMake's ability to set the correct
        # path of loki-transform.py if it was declared as an executable before (otherwise it
        # will assume it has been already on the path when CMake was called
        set( _LOKI_TRANSFORM loki-transform.py )
    endif()

    ####################################################
    ##  Idempotence mode:                             ##
    ##   * Internal "do-nothing" mode for Loki debug  ##
    ####################################################
    add_custom_command(
        OUTPUT loki-idem/cloudsc.idem.F90 loki-idem/cloudsc_driver_loki_mod.idem.F90
        COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/loki-idem
        COMMAND ${_LOKI_TRANSFORM} convert --mode idem --frontend ${LOKI_FRONTEND}
                      --config ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
                      --path ${CMAKE_CURRENT_SOURCE_DIR}
                      --header ${COMMON_MODULE}/yomphyder.F90
                      --cpp --include ${COMMON_INCLUDE}
                      --xmod ${XMOD_DIR}
                      --out-path ${CMAKE_CURRENT_BINARY_DIR}/loki-idem
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90
        COMMENT "[Loki] Pre-processing: mode=idem frontend=${LOKI_FRONTEND}"
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-idem
        SOURCES
            dwarf_cloudsc.F90
            loki-idem/cloudsc_driver_loki_mod.idem.F90
            loki-idem/cloudsc.idem.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries( dwarf-cloudsc-loki-idem PRIVATE cloudsc-common-lib )
    if( TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc-loki-idem PRIVATE OpenMP::OpenMP_Fortran )
    endif()


    ####################################################
    ##  SCA mode (Single Column Abstraction):         ##
    ##   * Extract de-vectorized SCA format code      ##
    ####################################################
    add_custom_command(
        OUTPUT loki-sca/cloudsc.sca.F90 loki-sca/cloudsc_driver_loki_mod.sca.F90
        COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/loki-sca
        COMMAND ${_LOKI_TRANSFORM} convert --mode sca --frontend ${LOKI_FRONTEND}
                      --config ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
                      --path ${CMAKE_CURRENT_SOURCE_DIR}
                      --header ${COMMON_MODULE}/yomphyder.F90
                      --cpp --include ${COMMON_INCLUDE}
                      --xmod ${XMOD_DIR}
                      --out-path ${CMAKE_CURRENT_BINARY_DIR}/loki-sca
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90
        COMMENT "[Loki] Pre-processing: mode=sca frontend=${LOKI_FRONTEND}"
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-sca
        SOURCES
            dwarf_cloudsc.F90
            loki-sca/cloudsc_driver_loki_mod.sca.F90
            loki-sca/cloudsc.sca.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries( dwarf-cloudsc-loki-sca PRIVATE cloudsc-common-lib )
    if( TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc-loki-sca PRIVATE OpenMP::OpenMP_Fortran )
    endif()


    ####################################################
    ##  CLAW-CPU mode:                                ##
    ##   * Generate SCA code with CLAW annotations    ##
    ##   * Process with CLAW (CPU layout and OpenMP)  ##
    ####################################################
    if( HAVE_CLOUDSC_LOKI_CLAW )
        add_custom_command(
            OUTPUT loki-claw-cpu/cloudsc.claw.cpu.F90 loki-claw-cpu/cloudsc_driver_loki_mod.claw.cpu.F90
            COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/loki-claw-cpu
            COMMAND ${_LOKI_TRANSFORM} convert --mode claw --frontend ${LOKI_FRONTEND}
                        --config ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
                        --path ${CMAKE_CURRENT_SOURCE_DIR}
                        --header ${COMMON_MODULE}/yomphyder.F90
                        --cpp --include ${COMMON_INCLUDE}
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

        ecbuild_add_executable( TARGET dwarf-cloudsc-loki-claw-cpu
            SOURCES
                dwarf_cloudsc.F90
                loki-claw-cpu/cloudsc_driver_loki_mod.claw.cpu.F90
                loki-claw-cpu/cloudsc.claw.cpu.F90
            DEFINITIONS ${CLOUDSC_DEFINITIONS}
        )
        target_link_libraries( dwarf-cloudsc-loki-claw-cpu PRIVATE cloudsc-common-lib )
        if( TARGET OpenMP::OpenMP_Fortran )
            target_link_libraries( dwarf-cloudsc-loki-claw-cpu PRIVATE OpenMP::OpenMP_Fortran )
        endif()
    endif()

    ####################################################
    ##  CLAW-GPU mode:                                ##
    ##   * Generate SCA code with CLAW annotations    ##
    ##   * Process with CLAW (GPU layout and OpenACC) ##
    ####################################################
    if( HAVE_CLOUDSC_LOKI_CLAW )
        add_custom_command(
            OUTPUT loki-claw-gpu/cloudsc.claw.gpu.F90 loki-claw-gpu/cloudsc_driver_loki_mod.claw.gpu.F90
            COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/loki-claw-gpu
            # Uses Loki-frontend CPP to switch to statement function variant again,
            # but suppresses inlining of stmt funcs by omitting `--include`
            COMMAND ${_LOKI_TRANSFORM} convert --mode claw --frontend ${LOKI_FRONTEND}
                          --config ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
                          --path ${CMAKE_CURRENT_SOURCE_DIR}
                          --header ${COMMON_MODULE}/yomphyder.F90
                          --cpp --include ${COMMON_INCLUDE}
                          --define CLOUDSC_GPU_TIMING ${CLOUDSC_DEFINE_STMT_FUNC}
                          --data-offload --remove-openmp
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

        ecbuild_add_executable( TARGET dwarf-cloudsc-loki-claw-gpu
            SOURCES
                dwarf_cloudsc.F90
                loki-claw-gpu/cloudsc_driver_loki_mod.claw.gpu.F90
                loki-claw-gpu/cloudsc.claw.gpu.F90
            DEFINITIONS ${CLOUDSC_DEFINITIONS}
        )
        if( TARGET OpenMP::OpenMP_Fortran )
            target_link_libraries( dwarf-cloudsc-loki-claw-gpu PRIVATE OpenMP::OpenMP_Fortran )
        endif()
        if( TARGET OpenACC::OpenACC_Fortran )
            target_link_libraries( dwarf-cloudsc-loki-claw-gpu PRIVATE OpenACC::OpenACC_Fortran )
            target_link_libraries( dwarf-cloudsc-loki-claw-gpu PRIVATE cloudsc-common-lib-static )
        else()
            target_link_libraries( dwarf-cloudsc-loki-claw-gpu PRIVATE cloudsc-common-lib )
        endif()
    endif()


    ####################################################
    ##  "Single Column Coalesced" (SCC) mode          ##
    ##   * Removes horizontal vector loops            ##
    ##   * Invokes compute kernel as `!$acc vector`   ##
    ####################################################
    add_custom_command(
        OUTPUT loki-scc/cloudsc.scc.F90 loki-scc/cloudsc_driver_loki_mod.scc.F90
        COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/loki-scc
        COMMAND ${_LOKI_TRANSFORM} convert --mode scc --frontend ${LOKI_FRONTEND}
                      --config ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
                      --path ${CMAKE_CURRENT_SOURCE_DIR}
                      --header ${COMMON_MODULE}/yomphyder.F90
                      --cpp --include ${COMMON_INCLUDE}
                      --define CLOUDSC_GPU_TIMING ${CLOUDSC_DEFINE_STMT_FUNC}
                      --data-offload --remove-openmp
                      --xmod ${XMOD_DIR}
                      --out-path ${CMAKE_CURRENT_BINARY_DIR}/loki-scc
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90
        COMMENT "[Loki] Pre-processing: mode=scc frontend=${LOKI_FRONTEND}"
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc
        SOURCES
            dwarf_cloudsc.F90
            loki-scc/cloudsc_driver_loki_mod.scc.F90
            loki-scc/cloudsc.scc.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    if( TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc-loki-scc PRIVATE OpenMP::OpenMP_Fortran )
    endif()
    if( TARGET OpenACC::OpenACC_Fortran )
        target_link_libraries( dwarf-cloudsc-loki-scc PRIVATE OpenACC::OpenACC_Fortran )
        target_link_libraries( dwarf-cloudsc-loki-scc PRIVATE cloudsc-common-lib-static )
    else()
        target_link_libraries( dwarf-cloudsc-loki-scc PRIVATE cloudsc-common-lib )
    endif()


    ####################################################
    ##  SCC-hoist mode                                ##
    ##   * SCC with vector loop hoisted               ##
    ##   * Kernel is "seq, but args are full blocks   ##
    ##   * Temporary arrays hoisted to driver         ##
    ####################################################
    add_custom_command(
        OUTPUT
            loki-scc-hoist/cloudsc.scc_hoist.F90
            loki-scc-hoist/cloudsc_driver_loki_mod.scc_hoist.F90
        COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-hoist
        COMMAND ${_LOKI_TRANSFORM} convert --mode scc-hoist --frontend ${LOKI_FRONTEND}
                      --config ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
                      --path ${CMAKE_CURRENT_SOURCE_DIR}
                      --header ${COMMON_MODULE}/yomphyder.F90
                      --cpp --include ${COMMON_INCLUDE}
                      --define CLOUDSC_GPU_TIMING ${CLOUDSC_DEFINE_STMT_FUNC}
                      --data-offload --remove-openmp
                      --xmod ${XMOD_DIR}
                      --out-path ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-hoist
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90
        COMMENT "[Loki] Pre-processing: mode=scc-hoist frontend=${LOKI_FRONTEND}"
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-hoist
        SOURCES
            dwarf_cloudsc.F90
            loki-scc-hoist/cloudsc_driver_loki_mod.scc_hoist.F90
            loki-scc-hoist/cloudsc.scc_hoist.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    if( TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc-loki-scc-hoist PRIVATE OpenMP::OpenMP_Fortran )
    endif()
    if( TARGET OpenACC::OpenACC_Fortran )
        target_link_libraries( dwarf-cloudsc-loki-scc-hoist PRIVATE OpenACC::OpenACC_Fortran )
        target_link_libraries( dwarf-cloudsc-loki-scc-hoist PRIVATE cloudsc-common-lib-static )
    else()
        target_link_libraries( dwarf-cloudsc-loki-scc-hoist PRIVATE cloudsc-common-lib )
    endif()


    ##############################################################################
    # C-transpilation mode for generating vectorized C host code (experimental!)
    ##############################################################################
    add_custom_command(
        OUTPUT loki-c/cloudsc_driver_loki_mod.c.F90
               loki-c/cloudsc_fc.F90 loki-c/cloudsc_c.c
               loki-c/yoethf_fc.F90 loki-c/yomcst_fc.F90
               loki-c/yoecldp_fc.F90
        COMMAND mkdir -p ${CMAKE_CURRENT_BINARY_DIR}/loki-c
        COMMAND ${_LOKI_TRANSFORM} transpile --frontend=${LOKI_FRONTEND}
                      --header ${COMMON_MODULE}/parkind1.F90
                      --header ${COMMON_MODULE}/yomphyder.F90
                      --header ${COMMON_MODULE}/yomcst.F90
                      --header ${COMMON_MODULE}/yoethf.F90
                      --header ${COMMON_MODULE}/yoecldp.F90
                      --header ${COMMON_MODULE}/fcttre_mod.F90
                      --header ${COMMON_MODULE}/fccld_mod.F90
                      --driver ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_driver_loki_mod.F90
                      --source ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc.F90
                      --cpp --include ${COMMON_INCLUDE}
                      --xmod ${XMOD_DIR}
                      --out-path ${CMAKE_CURRENT_BINARY_DIR}/loki-c
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90
        COMMENT "[Loki] Pre-processing: mode=transpile frontend=${LOKI_FRONTEND}"
    )

    # Define the CLAW-CPU build target for this variant
    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-c
        SOURCES
            dwarf_cloudsc.F90
            loki-c/cloudsc_driver_loki_mod.c.F90
            loki-c/cloudsc_fc.F90
            loki-c/cloudsc_c.c
            loki-c/yoethf_fc.F90
            loki-c/yomcst_fc.F90
            loki-c/yoecldp_fc.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries( dwarf-cloudsc-loki-c PRIVATE cloudsc-common-lib )
    if( TARGET OpenMP::OpenMP_Fortran AND TARGET OpenMP::OpenMP_C )
        target_link_libraries( dwarf-cloudsc-loki-c PRIVATE OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc-loki-c PRIVATE OpenMP::OpenMP_C )
    endif()

    # Create symlink for the input data
    if( HAVE_SERIALBOX )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )
    endif()

    if( HAVE_HDF5 )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../input.h5 )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../reference.h5 )
    endif()
endif()
