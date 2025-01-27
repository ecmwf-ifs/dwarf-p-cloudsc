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

if( HAVE_CLOUDSC_LOKI )

    ####################################################
    ##  Define various pre-processing modes via Loki  ##
    ####################################################

    set( COMMON_MODULE "${CMAKE_CURRENT_SOURCE_DIR}/../common/module" )
    set( COMMON_INCLUDE "${CMAKE_CURRENT_SOURCE_DIR}/../common/include" )
    set( XMOD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/xmod" )

    set( LOKI_FRONTEND "fp" CACHE STRING "Frontend parser for Loki transforms" )

    ####################################################
    ##  Idempotence mode:                             ##
    ##   * Internal "do-nothing" mode for Loki debug  ##
    ####################################################

    loki_transform(
        COMMAND convert
        OUTPUT
            loki-idem/cloudsc.idem.F90
            loki-idem/cloudsc_driver_loki_mod.idem.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-idem
        DEPENDS
            cloudsc.F90
            cloudsc_driver_loki_mod.F90
            ${_OMNI_DEPENDENCIES}
        MODE idem
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        CPP
        FRONTEND ${LOKI_FRONTEND}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        XMOD
            ${_TARGET_XMOD_DIR}
            ${XMOD_DIR}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-idem
        SOURCES
            dwarf_cloudsc.F90
            loki-idem/cloudsc_driver_loki_mod.idem.F90
            loki-idem/cloudsc.idem.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    # Set specific module directory to avoid aliasing of .mod files
    set_target_properties(  dwarf-cloudsc-loki-idem
        PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/loki-idem
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-serial
        COMMAND bin/dwarf-cloudsc-loki-idem
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-omp
        COMMAND bin/dwarf-cloudsc-loki-idem
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 4
        CONDITION HAVE_OMP
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-mpi
        COMMAND bin/dwarf-cloudsc-loki-idem
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        MPI 2
        OMP 1
        CONDITION HAVE_MPI
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-mpi-omp
        COMMAND bin/dwarf-cloudsc-loki-idem
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        MPI 2
        OMP 4
        CONDITION HAVE_OMP AND HAVE_MPI
    )

    ############################################################
    ##  Idempotence mode with pool-allocator for temporaries: ##
    ##   * Internal "do-nothing" mode for Loki debug          ##
    ############################################################

    loki_transform(
        COMMAND convert
        OUTPUT
            loki-idem-stack/cloudsc.idem_stack.F90
            loki-idem-stack/cloudsc_driver_loki_mod.idem_stack.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-idem-stack
        DEPENDS
            cloudsc.F90
            cloudsc_driver_loki_mod.F90
            ${_OMNI_DEPENDENCIES}
        MODE idem-stack
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        CPP
        FRONTEND ${LOKI_FRONTEND}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        XMOD
            ${_TARGET_XMOD_DIR}
            ${XMOD_DIR}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-idem-stack
        SOURCES
            dwarf_cloudsc.F90
            loki-idem-stack/cloudsc_driver_loki_mod.idem_stack.F90
            loki-idem-stack/cloudsc.idem_stack.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    # Set specific module directory to avoid aliasing of .mod files
    set_target_properties(  dwarf-cloudsc-loki-idem-stack
        PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/loki-idem-stack
    )


    if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        target_compile_options( dwarf-cloudsc-loki-idem-stack PRIVATE "-fcray-pointer" )
    elseif( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" OR CMAKE_Fortran_COMPILER_ID MATCHES "PGI" )
        target_compile_options( dwarf-cloudsc-loki-idem-stack PRIVATE "-Mcray=pointer" )
    endif()

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-stack-serial
        COMMAND bin/dwarf-cloudsc-loki-idem-stack
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-stack-omp
        COMMAND bin/dwarf-cloudsc-loki-idem-stack
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 4
        CONDITION HAVE_OMP
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-stack-mpi
        COMMAND bin/dwarf-cloudsc-loki-idem-stack
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        MPI 2
        OMP 1
        CONDITION HAVE_MPI
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-idem-stack-mpi-omp
        COMMAND bin/dwarf-cloudsc-loki-idem-stack
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        MPI 2
        OMP 4
        CONDITION HAVE_OMP AND HAVE_MPI
    )


    ####################################################
    ##  "Single Column Coalesced" (SCC) mode          ##
    ##   * Removes horizontal vector loops            ##
    ##   * Invokes compute kernel as `!$acc vector`   ##
    ####################################################

    loki_transform(
        COMMAND convert
        OUTPUT
            loki-scc/cloudsc.scc.F90
            loki-scc/cloudsc_driver_loki_mod.scc.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-scc
        DEPENDS
            cloudsc.F90
            cloudsc_driver_loki_mod.F90
            ${_OMNI_DEPENDENCIES}
        MODE scc
	CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
	# CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki_omp.config
        CPP
        DEFINITIONS
            CLOUDSC_GPU_TIMING
        FRONTEND ${LOKI_FRONTEND}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        XMOD
            ${_TARGET_XMOD_DIR}
            ${XMOD_DIR}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc
        SOURCES
            dwarf_cloudsc.F90
            loki-scc/cloudsc_driver_loki_mod.scc.F90
            loki-scc/cloudsc.scc.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    # Set specific module directory to avoid aliasing of .mod files
    set_target_properties(  dwarf-cloudsc-loki-scc
        PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/loki-scc
    )


    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-serial
        COMMAND bin/dwarf-cloudsc-loki-scc
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
        ENVIRONMENT "NVCOMPILER_ACC_CUDA_HEAPSIZE=128M"
    )

    ######################################################
    ##  "Single Column Coalesced" (SCC) mode with stack ##
    ##   * Removes horizontal vector loops              ##
    ##   * Invokes compute kernel as `!$acc vector`     ##
    ##   * Allocates temporaries using pool allocator   ##
    ######################################################

    loki_transform(
        COMMAND convert
        OUTPUT
            loki-scc-stack/cloudsc.scc_stack.F90
            loki-scc-stack/cloudsc_driver_loki_mod.scc_stack.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-stack
        DEPENDS
            cloudsc.F90
            cloudsc_driver_loki_mod.F90
            ${_OMNI_DEPENDENCIES}
        MODE scc-stack
	# CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
	CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki_omp.config
        CPP
        DEFINITIONS
            CLOUDSC_GPU_TIMING
        FRONTEND ${LOKI_FRONTEND}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        XMOD
            ${_TARGET_XMOD_DIR}
            ${XMOD_DIR}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-stack
        SOURCES
            dwarf_cloudsc.F90
            loki-scc-stack/cloudsc_driver_loki_mod.scc_stack.F90
            loki-scc-stack/cloudsc.scc_stack.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    # Set specific module directory to avoid aliasing of .mod files
    set_target_properties(  dwarf-cloudsc-loki-scc-stack
        PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-stack
    )


    if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
        target_compile_options( dwarf-cloudsc-loki-scc-stack PRIVATE "-fcray-pointer" )
    elseif( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" OR CMAKE_Fortran_COMPILER_ID MATCHES "PGI" )
        target_compile_options( dwarf-cloudsc-loki-scc-stack PRIVATE "-Mcray=pointer" )
    endif()

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-stack-serial
        COMMAND bin/dwarf-cloudsc-loki-scc-stack
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )


    ####################################################
    ##  SCC-hoist mode                                ##
    ##   * SCC with vector loop hoisted               ##
    ##   * Kernel is "seq, but args are full blocks   ##
    ##   * Temporary arrays hoisted to driver         ##
    ####################################################

    loki_transform(
        COMMAND convert
        OUTPUT
            loki-scc-hoist/cloudsc.scc_hoist.F90
            loki-scc-hoist/cloudsc_driver_loki_mod.scc_hoist.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-hoist
        DEPENDS
            cloudsc.F90
            cloudsc_driver_loki_mod.F90
            ${_OMNI_DEPENDENCIES}
        MODE scc-hoist
	CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
	# CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki_omp.config
        CPP
        DEFINITIONS
            CLOUDSC_GPU_TIMING
        FRONTEND ${LOKI_FRONTEND}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        XMOD
            ${_TARGET_XMOD_DIR}
            ${XMOD_DIR}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-hoist
        SOURCES
            dwarf_cloudsc.F90
            loki-scc-hoist/cloudsc_driver_loki_mod.scc_hoist.F90
            loki-scc-hoist/cloudsc.scc_hoist.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    # Set specific module directory to avoid aliasing of .mod files
    set_target_properties(  dwarf-cloudsc-loki-scc-hoist
        PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-hoist
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-hoist-serial
        COMMAND bin/dwarf-cloudsc-loki-scc-hoist
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )


    ####################################################
    ##  "Single Column Coalesced Field" (SCC-FIELD)   ##
    ##   * Removes horizontal vector loops            ##
    ##   * Invokes compute kernel as `!$acc vector`   ##
    ##   * Uses Field API for device data offloading  ##
    ####################################################
if ( HAVE_FIELD_API AND field_api_HAVE_ACC )
    loki_transform(
        COMMAND convert
        OUTPUT
            loki-scc-field/cloudsc.scc_field.F90
            loki-scc-field/cloudsc_driver_field_loki_mod.scc_field.F90
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-field
        DEPENDS
            cloudsc.F90
            cloudsc_driver_field_loki_mod.F90
            ${_OMNI_DEPENDENCIES}
        MODE scc-field
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki_field_offload.config
        CPP
        DEFINITIONS
            CLOUDSC_GPU_TIMING
            ${CLOUDSC_DEFINE_STMT_FUNC}
        FRONTEND ${LOKI_FRONTEND}
        HEADERS
            ${COMMON_MODULE}/yomcst.F90
            ${COMMON_MODULE}/yomphyder.F90
            ${COMMON_MODULE}/yoethf.F90
            ${COMMON_MODULE}/yoecldp.F90
            ${COMMON_MODULE}/cloudsc_field_state_mod.F90
            ${COMMON_MODULE}/cloudsc_flux_type_mod.F90
            ${COMMON_MODULE}/cloudsc_aux_type_mod.F90
            ${COMMON_MODULE}/cloudsc_state_type_mod.F90
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        XMOD
            ${_TARGET_XMOD_DIR}
            ${XMOD_DIR}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-field
        SOURCES
            dwarf_cloudsc.F90
            loki-scc-field/cloudsc.scc_field.F90
            loki-scc-field/cloudsc_driver_field_loki_mod.scc_field.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_FIELD
    )
    # Set specific module directory to avoid aliasing of .mod files
    set_target_properties(  dwarf-cloudsc-loki-scc-field
    	PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-field
    )


    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-field
        COMMAND bin/dwarf-cloudsc-loki-scc-field
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
        ENVIRONMENT "NVCOMPILER_ACC_CUDA_HEAPSIZE=128M"
    )

endif()


    ####################################################
    ##  SCC CUF (CUDA Fortran)                        ##
    ##   * SCC with CUDA Fortran (CUF)                ##
    ##   * allows for 3 types/kind of transformation  ##
    ##    * parametrised (dimensions) local arrays    ##
    ##    * host side hoisted local arrays            ##
    ##    * dynamic memory allocation for local       ##
    ##      arrays (not recommended)                  ##
    ####################################################


if( HAVE_CUDA )

    # scc-cuf-parametrise
    loki_transform(
        COMMAND convert
        MODE cuf-parametrise
        FRONTEND ${LOKI_FRONTEND}
        CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
        SOURCES ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS
            ${COMMON_MODULE}/yomcst.F90
            ${COMMON_MODULE}/yomphyder.F90
            ${COMMON_MODULE}/yoethf.F90
            ${COMMON_MODULE}/yoecldp.F90
        INCLUDES ${COMMON_INCLUDE}
        DEFINITIONS CLOUDSC_GPU_TIMING
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-cuf-parametrise
        OUTPUT
            loki-scc-cuf-parametrise/cuf_cloudsc_driver_loki_mod.cuf_parametrise.F90
            loki-scc-cuf-parametrise/cuf_cloudsc.cuf_parametrise.F90
        DEPENDS cuf_cloudsc.F90 cuf_cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )

    set_source_files_properties(
        loki-scc-cuf-parametrise/cuf_cloudsc_driver_loki_mod.cuf_parametrise.F90
        loki-scc-cuf-parametrise/cuf_cloudsc.cuf_parametrise.F90
        PROPERTIES COMPILE_FLAGS "-cuda -gpu=maxregcount:128"
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-cuf-parametrise
        SOURCES
            dwarf_cloudsc.F90
            loki-scc-cuf-parametrise/cuf_cloudsc_driver_loki_mod.cuf_parametrise.F90
            loki-scc-cuf-parametrise/cuf_cloudsc.cuf_parametrise.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_CUF
    )
    # Set specific module directory to avoid aliasing of .mod files
    set_target_properties(  dwarf-cloudsc-loki-scc-cuf-parametrise
        PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-cuf-parametrise
    )

    target_link_options(dwarf-cloudsc-loki-scc-cuf-parametrise PUBLIC "-cuda")

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-cuf-parametrise-serial
        COMMAND bin/dwarf-cloudsc-loki-scc-cuf-parametrise
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )

    # scc-cuf-hoist
    loki_transform(
        COMMAND convert
        MODE cuf-hoist
        FRONTEND ${LOKI_FRONTEND}
        CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
        SOURCES ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS
            ${COMMON_MODULE}/yomcst.F90
            ${COMMON_MODULE}/yomphyder.F90
            ${COMMON_MODULE}/yoethf.F90
            ${COMMON_MODULE}/yoecldp.F90
        INCLUDES ${COMMON_INCLUDE}
        DEFINITIONS CLOUDSC_GPU_TIMING
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-cuf-hoist
        OUTPUT
            loki-scc-cuf-hoist/cuf_cloudsc_driver_loki_mod.cuf_hoist.F90
            loki-scc-cuf-hoist/cuf_cloudsc.cuf_hoist.F90
        DEPENDS cuf_cloudsc.F90 cuf_cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )

    set_source_files_properties(
        loki-scc-cuf-hoist/cuf_cloudsc_driver_loki_mod.cuf_hoist.F90
        loki-scc-cuf-hoist/cuf_cloudsc.cuf_hoist.F90
        PROPERTIES COMPILE_FLAGS "-cuda -gpu=maxregcount:128"
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-cuf-hoist
        SOURCES
            dwarf_cloudsc.F90
            loki-scc-cuf-hoist/cuf_cloudsc_driver_loki_mod.cuf_hoist.F90
            loki-scc-cuf-hoist/cuf_cloudsc.cuf_hoist.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_CUF
    )
    # Set specific module directory to avoid aliasing of .mod files
    set_target_properties(  dwarf-cloudsc-loki-scc-cuf-hoist
        PROPERTIES Fortran_MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-cuf-hoist
    )

    target_link_options(dwarf-cloudsc-loki-scc-cuf-hoist PUBLIC "-cuda")

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-cuf-hoist-serial
        COMMAND bin/dwarf-cloudsc-loki-scc-cuf-hoist
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )

endif()

    #############################################################################
    ## C-transpilation mode for generating vectorized C host code (experimental!)
    #############################################################################

    if ( NOT HAVE_SINGLE_PRECISION )

        loki_transform(
            COMMAND convert
            MODE c FRONTEND ${LOKI_FRONTEND} CPP
            CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
            SOURCES ${CMAKE_CURRENT_SOURCE_DIR}
            HEADERS
                ${COMMON_MODULE}/parkind1.F90
                ${COMMON_MODULE}/yomphyder.F90
                ${COMMON_MODULE}/yomcst.F90
                ${COMMON_MODULE}/yoethf.F90
                ${COMMON_MODULE}/yoecldp.F90
            INCLUDES ${COMMON_INCLUDE}
            XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
            BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-c
            OUTPUT
                loki-c/cloudsc_driver_loki_mod.c.F90
                loki-c/cloudsc_fc.F90 loki-c/cloudsc_c.c
                loki-c/yoethf_fc.F90 loki-c/yomcst_fc.F90
                loki-c/yoecldp_fc.F90
            DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
        )

        # Define the build target for the Loki-C variant
        ecbuild_add_executable( TARGET dwarf-cloudsc-loki-c
            SOURCES
                dwarf_cloudsc.F90
                loki-c/cloudsc_driver_loki_mod.c.F90
                loki-c/cloudsc_fc.F90
                loki-c/cloudsc_c.c
                loki-c/yoethf_fc.F90
                loki-c/yomcst_fc.F90
                loki-c/yoecldp_fc.F90
            LIBS
                cloudsc-common-lib
            DEFINITIONS ${CLOUDSC_DEFINITIONS}
        )

        ecbuild_add_test(
            TARGET dwarf-cloudsc-loki-c-serial
            COMMAND bin/dwarf-cloudsc-loki-c
            ARGS 1 100 16
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
            OMP 1
        )
        ecbuild_add_test(
            TARGET dwarf-cloudsc-loki-c-omp
            COMMAND bin/dwarf-cloudsc-loki-c
            ARGS 4 100 16
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
            OMP 4
            CONDITION HAVE_OMP
        )

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
