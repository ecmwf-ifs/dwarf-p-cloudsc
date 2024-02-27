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

# Define the CLAW-based Loki dwarf variants as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_LOKI_CLAW
    DESCRIPTION "Build optimized CPU and GPU version of CLOUDSC derived from Loki+CLAW" DEFAULT OFF
    CONDITION HAVE_CLOUDSC_LOKI
)

function( cloudsc_xmod _TARGET )

    if( TARGET clawfc AND ${LOKI_FRONTEND} STREQUAL "omni" )

        # Ugly hack: OMNI needs the xmod-file for cloudsc.F90 to be able to
        # parse the driver file successfully. However, the scheduler currently
        # doesn't take this into account and fails when parsing driver before
        # kernel file.
        # (Note: the problem vanishes in serial builds as there the C-transpile
        # target is built first which doesn't use the scheduler and therefore
        # creates the necessary xmod files for us)
        # TODO: This can be removed once the scheduler is aware of these dependencies
        # and parses files in the right order

        set( _TARGET_XMOD_DIR "${CMAKE_CURRENT_BINARY_DIR}/${_TARGET}" )
        set( _TARGET_XMOD_DIR ${_TARGET_XMOD_DIR} PARENT_SCOPE )
        file( MAKE_DIRECTORY ${_TARGET_XMOD_DIR} )

        generate_xmod(
            OUTPUT ${_TARGET_XMOD_DIR}/cloudsc.xml
            SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc.F90
            XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        )

        # Use XML files as dependencies (not xmod) as they are updated by later calls of
        # F_Front (and thus would trigger new execution rounds)
        set( _OMNI_DEPENDENCIES ${_TARGET_XMOD_DIR}/cloudsc.xml PARENT_SCOPE )

    else()

        set( _TARGET_XMOD_DIR "" PARENT_SCOPE)
        set( _OMNI_DEPENDENCIES "" PARENT_SCOPE )

    endif()

endfunction()

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
      set( CLOUDSC_DEFINE_STMT_FUNC CLOUDSC_STMT_FUNC )
    endif()

    ####################################################
    ##  Idempotence mode:                             ##
    ##   * Internal "do-nothing" mode for Loki debug  ##
    ####################################################

 if (FALSE)

    cloudsc_xmod( loki-idem )

    loki_transform_convert(
        MODE idem FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS ${COMMON_MODULE}/yomphyder.F90
        INCLUDES ${COMMON_INCLUDE}
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-idem
        OUTPUT loki-idem/cloudsc.idem.F90 loki-idem/cloudsc_driver_loki_mod.idem.F90
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
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

    cloudsc_xmod( loki-idem-stack )

    loki_transform_convert(
        MODE idem-stack FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS ${COMMON_MODULE}/yomphyder.F90
        INCLUDES ${COMMON_INCLUDE}
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-idem-stack
        OUTPUT loki-idem-stack/cloudsc.idem_stack.F90 loki-idem-stack/cloudsc_driver_loki_mod.idem_stack.F90
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-idem-stack
        SOURCES
            dwarf_cloudsc.F90
            stack_mod.F90
            loki-idem-stack/cloudsc_driver_loki_mod.idem_stack.F90
            loki-idem-stack/cloudsc.idem_stack.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
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
    ##  SCA mode (Single Column Abstraction):         ##
    ##   * Extract de-vectorized SCA format code      ##
    ####################################################

    cloudsc_xmod( loki-sca )

    loki_transform_convert(
        MODE sca FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS ${COMMON_MODULE}/yomphyder.F90
        INCLUDES ${COMMON_INCLUDE}
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-sca
        OUTPUT loki-sca/cloudsc.sca.F90 loki-sca/cloudsc_driver_loki_mod.sca.F90
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-sca
        SOURCES
            dwarf_cloudsc.F90
            loki-sca/cloudsc_driver_loki_mod.sca.F90
            loki-sca/cloudsc.sca.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-sca-serial
        COMMAND bin/dwarf-cloudsc-loki-sca
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-sca-mpi
        COMMAND bin/dwarf-cloudsc-loki-sca
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        MPI 2
        OMP 1
        CONDITION HAVE_MPI
    )

    ####################################################
    ##  CLAW-CPU mode:                                ##
    ##   * Generate SCA code with CLAW annotations    ##
    ##   * Process with CLAW (CPU layout and OpenMP)  ##
    ####################################################
    if( HAVE_CLOUDSC_LOKI_CLAW AND TARGET clawfc )

        cloudsc_xmod( loki-claw-cpu )

        loki_transform_convert(
            MODE claw FRONTEND ${LOKI_FRONTEND} CPP
            CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
            PATH ${CMAKE_CURRENT_SOURCE_DIR}
            HEADERS ${COMMON_MODULE}/yomphyder.F90
            INCLUDES ${COMMON_INCLUDE}
            XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
            OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-claw-cpu
            OUTPUT loki-claw-cpu/cloudsc.claw.F90 loki-claw-cpu/cloudsc_driver_loki_mod.claw.F90
            DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
        )

        # We purposefully suppress CLAWs insertion of OpenMP loops,
        # as they are already present in the outer driver.
        claw_compile(
            MODEL_CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/claw_cloudsc.config
            TARGET cpu DIRECTIVE none
            INCLUDES ${COMMON_INCLUDE}
            XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
            SOURCE loki-claw-cpu/cloudsc.claw.F90
            OUTPUT loki-claw-cpu/cloudsc.claw.cpu.F90
        )
        claw_compile(
            MODEL_CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/claw_cloudsc.config
            TARGET cpu DIRECTIVE none
            INCLUDES ${COMMON_INCLUDE}
            XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
            SOURCE loki-claw-cpu/cloudsc_driver_loki_mod.claw.F90
            OUTPUT loki-claw-cpu/cloudsc_driver_loki_mod.claw.cpu.F90
            DEPENDS loki-claw-cpu/cloudsc.claw.cpu.F90
        )

        ecbuild_add_executable( TARGET dwarf-cloudsc-loki-claw-cpu
            SOURCES
                dwarf_cloudsc.F90
                loki-claw-cpu/cloudsc_driver_loki_mod.claw.cpu.F90
                loki-claw-cpu/cloudsc.claw.cpu.F90
            LIBS
                cloudsc-common-lib
            DEFINITIONS ${CLOUDSC_DEFINITIONS}
        )

        ecbuild_add_test(
            TARGET dwarf-cloudsc-loki-claw-cpu-serial
            COMMAND bin/dwarf-cloudsc-loki-claw-cpu
            ARGS 1 100 16
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
            OMP 1
        )
        ecbuild_add_test(
            TARGET dwarf-cloudsc-loki-claw-cpu-omp
            COMMAND bin/dwarf-cloudsc-loki-claw-cpu
            ARGS 4 100 16
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
            OMP 4
            CONDITION HAVE_OMP
        )
        ecbuild_add_test(
            TARGET dwarf-cloudsc-loki-claw-cpu-mpi
            COMMAND bin/dwarf-cloudsc-loki-claw-cpu
            ARGS 1 100 16
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
            MPI 2
            OMP 1
            CONDITION HAVE_MPI
        )
        ecbuild_add_test(
            TARGET dwarf-cloudsc-loki-claw-cpu-mpi-omp
            COMMAND bin/dwarf-cloudsc-loki-claw-cpu
            ARGS 4 100 16
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
            MPI 2
            OMP 4
            CONDITION HAVE_OMP AND HAVE_MPI
        )
    endif()

    ####################################################
    ##  CLAW-GPU mode:                                ##
    ##   * Generate SCA code with CLAW annotations    ##
    ##   * Process with CLAW (GPU layout and OpenACC) ##
    ####################################################
    if( HAVE_CLOUDSC_LOKI_CLAW AND TARGET clawfc )

        cloudsc_xmod( loki-claw-gpu )

        # Uses Loki-frontend CPP to switch to statement function variant again,
        # but suppresses inlining of stmt funcs by omitting `--include`
        loki_transform_convert(
            MODE claw FRONTEND ${LOKI_FRONTEND} CPP
            CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
            PATH ${CMAKE_CURRENT_SOURCE_DIR}
            HEADERS ${COMMON_MODULE}/yomphyder.F90
            INCLUDES ${COMMON_INCLUDE}
            DEFINITIONS CLOUDSC_GPU_TIMING ${CLOUDSC_DEFINE_STMT_FUNC}
            DATA_OFFLOAD REMOVE_OPENMP
            XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
            OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-claw-gpu
            OUTPUT loki-claw-gpu/cloudsc.claw.F90 loki-claw-gpu/cloudsc_driver_loki_mod.claw.F90
            DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
        )

        claw_compile(
            MODEL_CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/claw_cloudsc.config
            TARGET gpu DIRECTIVE openacc
            INCLUDES ${COMMON_INCLUDE}
            XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
            SOURCE loki-claw-gpu/cloudsc.claw.F90
            OUTPUT loki-claw-gpu/cloudsc.claw.gpu.F90
        )
        claw_compile(
            MODEL_CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/claw_cloudsc.config
            TARGET gpu DIRECTIVE openacc
            INCLUDES ${COMMON_INCLUDE}
            XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
            SOURCE loki-claw-gpu/cloudsc_driver_loki_mod.claw.F90
            OUTPUT loki-claw-gpu/cloudsc_driver_loki_mod.claw.gpu.F90
            DEPENDS loki-claw-gpu/cloudsc.claw.gpu.F90
        )

        ecbuild_add_executable( TARGET dwarf-cloudsc-loki-claw-gpu
            SOURCES
                dwarf_cloudsc.F90
                loki-claw-gpu/cloudsc_driver_loki_mod.claw.gpu.F90
                loki-claw-gpu/cloudsc.claw.gpu.F90
            LIBS
                cloudsc-common-lib
            DEFINITIONS ${CLOUDSC_DEFINITIONS}
        )

        ecbuild_add_test(
            TARGET dwarf-cloudsc-loki-claw-gpu-serial
            COMMAND bin/dwarf-cloudsc-loki-claw-gpu
            ARGS 1 1280 128
            WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
            OMP 1
        )
    endif()


endif() # FALSE

if (FALSE)
    ####################################################
    ##  "Single Column Coalesced" (SCC) mode          ##
    ##   * Removes horizontal vector loops            ##
    ##   * Invokes compute kernel as `!$acc vector`   ##
    ####################################################

    cloudsc_xmod( loki-scc )

    loki_transform_convert(
        MODE scc FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS ${COMMON_MODULE}/yomphyder.F90 ${COMMON_MODULE}/yoecldp.F90
        INCLUDES ${COMMON_INCLUDE}
        DEFINITIONS CLOUDSC_GPU_TIMING ${CLOUDSC_DEFINE_STMT_FUNC}
        DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-scc
        OUTPUT loki-scc/cloudsc.scc.F90 loki-scc/cloudsc_driver_loki_mod.scc.F90
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
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

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-serial
        COMMAND bin/dwarf-cloudsc-loki-scc
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
        ENVIRONMENT "NVCOMPILER_ACC_CUDA_HEAPSIZE=128M"
    )

# if (FALSE)

    # scc-cpu
    cloudsc_xmod( loki-scc-cpu )

    loki_transform_convert(
        MODE scc-cpu FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS
            ${COMMON_MODULE}/yomcst.F90
            ${COMMON_MODULE}/yomphyder.F90
            ${COMMON_MODULE}/yoethf.F90
            ${COMMON_MODULE}/yoecldp.F90
        INCLUDES ${COMMON_INCLUDE}
	# DEFINITIONS CLOUDSC_GPU_TIMING
        DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-cpu
        OUTPUT
            loki-scc-cpu/cuf_cloudsc_driver_loki_mod.scc_cpu.F90
            loki-scc-cpu/cuf_cloudsc.scc_cpu.F90
            DEPENDS cuf_cloudsc.F90 cuf_cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-cpu
        SOURCES
            dwarf_cloudsc.F90
            loki-scc-cpu/cuf_cloudsc_driver_loki_mod.scc_cpu.F90
            loki-scc-cpu/cuf_cloudsc.scc_cpu.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_CUF
    )

endif()
##


if( HAVE_CUDA )
if (TRUE)
    # scc-cuf-parametrise
    cloudsc_xmod( loki-scc-cuf-parametrise-new )

    loki_transform_convert(
        MODE cuf-parametrise-new FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS
            ${COMMON_MODULE}/yomcst.F90
            ${COMMON_MODULE}/yomphyder.F90
            ${COMMON_MODULE}/yoethf.F90
            ${COMMON_MODULE}/yoecldp.F90
	    ${COMMON_MODULE}/fcttre_mod.F90
            ${COMMON_MODULE}/fccld_mod.F90
        INCLUDES ${COMMON_INCLUDE}
	DEFINITIONS CLOUDSC_GPU_TIMING CUF_IMPORTS 
        DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-cuf-parametrise-new
        OUTPUT
            loki-scc-cuf-parametrise-new/cuf_cloudsc_driver_loki_mod.cuf_parametrise_new.F90
            loki-scc-cuf-parametrise-new/cuf_cloudsc.cuf_parametrise_new.F90
            DEPENDS cuf_cloudsc.F90 cuf_cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )

    set_source_files_properties(
        loki-scc-cuf-parametrise-new/cuf_cloudsc_driver_loki_mod.cuf_parametrise_new.F90
        loki-scc-cuf-parametrise-new/cuf_cloudsc.cuf_parametrise_new.F90
	# PROPERTIES COMPILE_FLAGS "-g -cuda -gpu=debug,maxregcount:128 -Minfo"
	PROPERTIES COMPILE_FLAGS "-cuda -gpu=maxregcount:128 -Minfo"
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-cuf-parametrise-new
        SOURCES
            dwarf_cloudsc.F90
            loki-scc-cuf-parametrise-new/cuf_cloudsc_driver_loki_mod.cuf_parametrise_new.F90
            loki-scc-cuf-parametrise-new/cuf_cloudsc.cuf_parametrise_new.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_CUF
    )

    target_link_options(dwarf-cloudsc-loki-scc-cuf-parametrise-new PUBLIC "-cuda")

    # scc-cuf-hoist
    cloudsc_xmod( loki-scc-cuf-hoist-new )

    loki_transform_convert(
        MODE cuf-hoist-new FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS
            ${COMMON_MODULE}/yomcst.F90
            ${COMMON_MODULE}/yomphyder.F90
            ${COMMON_MODULE}/yoethf.F90
            ${COMMON_MODULE}/yoecldp.F90
	    ${COMMON_MODULE}/fcttre_mod.F90
            ${COMMON_MODULE}/fccld_mod.F90
        INCLUDES ${COMMON_INCLUDE}
	DEFINITIONS CLOUDSC_GPU_TIMING CUF_IMPORTS
        DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-cuf-hoist-new
        OUTPUT
            loki-scc-cuf-hoist-new/cuf_cloudsc_driver_loki_mod.cuf_hoist_new.F90
            loki-scc-cuf-hoist-new/cuf_cloudsc.cuf_hoist_new.F90
            DEPENDS cuf_cloudsc.F90 cuf_cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )

    set_source_files_properties(
        loki-scc-cuf-hoist-new/cuf_cloudsc_driver_loki_mod.cuf_hoist_new.F90
        loki-scc-cuf-hoist-new/cuf_cloudsc.cuf_hoist_new.F90
	# PROPERTIES COMPILE_FLAGS "-g -cuda -gpu=debug,maxregcount:128 -Minfo"
	PROPERTIES COMPILE_FLAGS "-cuda -gpu=maxregcount:128 -Minfo"
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-cuf-hoist-new
        SOURCES
            dwarf_cloudsc.F90
            loki-scc-cuf-hoist-new/cuf_cloudsc_driver_loki_mod.cuf_hoist_new.F90
            loki-scc-cuf-hoist-new/cuf_cloudsc.cuf_hoist_new.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_CUF
    )

    target_link_options(dwarf-cloudsc-loki-scc-cuf-hoist-new PUBLIC "-cuda")
endif()
###

if (TRUE)
    enable_language(CUDA)
    set(CMAKE_CUDA_FLAGS "-O3 -ffast-math")
    find_package(CUDAToolkit)
    link_directories("/usr/local/apps/nvidia/22.11/Linux_x86_64/22.11/math_libs/lib64")

    #### CUDA PARAMETRISE #####
    cloudsc_xmod( loki-c-scc-parametrise )

    loki_transform_convert(
        MODE c-parametrise FRONTEND ${LOKI_FRONTEND} CPP
	CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
	# CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
	PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS
                ${COMMON_MODULE}/parkind1.F90
                ${COMMON_MODULE}/yomphyder.F90
                ${COMMON_MODULE}/yomcst.F90
                ${COMMON_MODULE}/yoethf.F90
                ${COMMON_MODULE}/yoecldp.F90
                ${COMMON_MODULE}/fcttre_mod.F90
                ${COMMON_MODULE}/fccld_mod.F90
        INCLUDES ${COMMON_INCLUDE}
	# DEFINITIONS CLOUDSC_GPU_TIMING ${CLOUDSC_DEFINE_STMT_FUNC}
        DEFINITIONS CLOUDSC_GPU_TIMING
        DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
	OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-c-scc-parametrise
        OUTPUT
            loki-c-scc-parametrise/cuf_cloudsc_driver_loki_mod.c_parametrise.F90
            loki-c-scc-parametrise/cloudsc_cuf_fc.F90 loki-c-scc-parametrise/cloudsc_cuf_c.c
            loki-c-scc-parametrise/yoethf_fc.F90 loki-c-scc-parametrise/yomcst_fc.F90
            loki-c-scc-parametrise/yoecldp_fc.F90
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )


      ecbuild_add_library(
        TARGET dwarf-cloudsc-scc-parametrise-gpu-lib
        INSTALL_HEADERS LISTED
        SOURCES
	    # yoecldp_c.h
            loki-c-scc-parametrise/cloudsc_cuf_c.c
    	    loki-c-scc-parametrise/cloudsc_cuf_fc.F90
	    PUBLIC_INCLUDES
    	    	$<INSTALL_INTERFACE:include>
            	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
    	    PUBLIC_LIBS
    	        cloudsc-common-lib
    )

    target_include_directories(dwarf-cloudsc-scc-parametrise-gpu-lib PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

    set_source_files_properties(loki-c-scc-parametrise/cloudsc_cuf_c.c PROPERTIES LANGUAGE CUDA)

   target_compile_options(dwarf-cloudsc-scc-parametrise-gpu-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
	   --generate-line-info
	   >)


      cloudsc_add_compile_options(
	 SOURCES dwarf_cloudsc.F90
         FLAGS   "-cuda") # -gpu=maxrregcount ...
    target_include_directories(
        dwarf-cloudsc-scc-parametrise-gpu-lib
        PUBLIC
        ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    )

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-loki-scc-parametrise-cuda
	SOURCES dwarf_cloudsc.F90 loki-c-scc-parametrise/cuf_cloudsc_driver_loki_mod.c_parametrise.F90
	LIBS dwarf-cloudsc-scc-parametrise-gpu-lib
    )

    target_compile_definitions(dwarf-cloudsc-loki-scc-parametrise-cuda PUBLIC CUF_DRIVER=1)

    set_target_properties(dwarf-cloudsc-loki-scc-parametrise-cuda PROPERTIES LINKER_LANGUAGE Fortran)

    set_target_properties( dwarf-cloudsc-loki-scc-parametrise-cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_link_libraries(dwarf-cloudsc-loki-scc-parametrise-cuda dwarf-cloudsc-scc-parametrise-gpu-lib cudadevrt cudart)

    #### CUDA HOIST #####
    cloudsc_xmod( loki-c-scc-hoist )

    loki_transform_convert(
        MODE c-hoist FRONTEND ${LOKI_FRONTEND} CPP
	CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
	# CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
	PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS
                ${COMMON_MODULE}/parkind1.F90
                ${COMMON_MODULE}/yomphyder.F90
                ${COMMON_MODULE}/yomcst.F90
                ${COMMON_MODULE}/yoethf.F90
                ${COMMON_MODULE}/yoecldp.F90
                ${COMMON_MODULE}/fcttre_mod.F90
                ${COMMON_MODULE}/fccld_mod.F90
        INCLUDES ${COMMON_INCLUDE}
	# DEFINITIONS CLOUDSC_GPU_TIMING ${CLOUDSC_DEFINE_STMT_FUNC}
        DEFINITIONS CLOUDSC_GPU_TIMING
        DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
	OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-c-scc-hoist
        OUTPUT
            loki-c-scc-hoist/cuf_cloudsc_driver_loki_mod.c_hoist.F90
            loki-c-scc-hoist/cloudsc_cuf_fc.F90 loki-c-scc-hoist/cloudsc_cuf_c.c
            loki-c-scc-hoist/yoethf_fc.F90 loki-c-scc-hoist/yomcst_fc.F90
            loki-c-scc-hoist/yoecldp_fc.F90
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )


      ecbuild_add_library(
        TARGET dwarf-cloudsc-scc-hoist-gpu-lib
        INSTALL_HEADERS LISTED
        SOURCES
	    # yoecldp_c.h
            loki-c-scc-hoist/cloudsc_cuf_c.c
    	    loki-c-scc-hoist/cloudsc_cuf_fc.F90
	    PUBLIC_INCLUDES
    	    	$<INSTALL_INTERFACE:include>
            	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
    	    PUBLIC_LIBS
    	        cloudsc-common-lib
    )

    target_include_directories(dwarf-cloudsc-scc-hoist-gpu-lib PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

    set_source_files_properties(loki-c-scc-hoist/cloudsc_cuf_c.c PROPERTIES LANGUAGE CUDA)

   target_compile_options(dwarf-cloudsc-scc-hoist-gpu-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
	   --generate-line-info
	   >)


      cloudsc_add_compile_options(
	 SOURCES dwarf_cloudsc.F90
         FLAGS   "-cuda") # -gpu=maxrregcount ...
    target_include_directories(
        dwarf-cloudsc-scc-hoist-gpu-lib
        PUBLIC
        ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    )

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-loki-scc-hoist-cuda
	SOURCES dwarf_cloudsc.F90 loki-c-scc-hoist/cuf_cloudsc_driver_loki_mod.c_hoist.F90
	LIBS dwarf-cloudsc-scc-hoist-gpu-lib
    )

    target_compile_definitions(dwarf-cloudsc-loki-scc-hoist-cuda PUBLIC CUF_DRIVER=1)

    set_target_properties(dwarf-cloudsc-loki-scc-hoist-cuda PROPERTIES LINKER_LANGUAGE Fortran)

    set_target_properties( dwarf-cloudsc-loki-scc-hoist-cuda PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_link_libraries(dwarf-cloudsc-loki-scc-hoist-cuda dwarf-cloudsc-scc-hoist-gpu-lib cudadevrt cudart)

endif()

###
endif()
##

if (HAVE_HIP)

    enable_language(HIP)
    set(CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} -O3 -ffast-math")

    set(CMAKE_C_COMPILER "${ROCM_PATH}/bin/hipcc")
    set(CMAKE_CXX_COMPILER "${ROCM_PATH}/bin/hipcc")

    # PARAMETRISE
    cloudsc_xmod( loki-hip-scc-parametrise )

    loki_transform_convert(
        MODE hip-parametrise FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
        # CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS
                ${COMMON_MODULE}/parkind1.F90
                ${COMMON_MODULE}/yomphyder.F90
                ${COMMON_MODULE}/yomcst.F90
                ${COMMON_MODULE}/yoethf.F90
                ${COMMON_MODULE}/yoecldp.F90
                ${COMMON_MODULE}/fcttre_mod.F90
                ${COMMON_MODULE}/fccld_mod.F90
        INCLUDES ${COMMON_INCLUDE}
        # DEFINITIONS CLOUDSC_GPU_TIMING ${CLOUDSC_DEFINE_STMT_FUNC}
        DEFINITIONS CLOUDSC_GPU_TIMING
        DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-hip-scc-parametrise
        OUTPUT
            loki-hip-scc-parametrise/cuf_cloudsc_driver_loki_mod.hip_parametrise.F90
            loki-hip-scc-parametrise/cloudsc_cuf_fc.F90 loki-hip-scc-parametrise/cloudsc_cuf_c.c
            loki-hip-scc-parametrise/yoethf_fc.F90 loki-hip-scc-parametrise/yomcst_fc.F90
            loki-hip-scc-parametrise/yoecldp_fc.F90
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )


      ecbuild_add_library(
        TARGET dwarf-cloudsc-scc-parametrise-gpu-lib
        INSTALL_HEADERS LISTED
        SOURCES
            # yoecldp_c.h
            loki-hip-scc-parametrise/cloudsc_cuf_c.c
            loki-hip-scc-parametrise/cloudsc_cuf_fc.F90
            PUBLIC_INCLUDES
                $<INSTALL_INTERFACE:include>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
            PUBLIC_LIBS
                cloudsc-common-lib
    )

   target_include_directories(dwarf-cloudsc-scc-parametrise-gpu-lib PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

   set_source_files_properties(loki-hip-scc-parametrise/cloudsc_cuf_c.c PROPERTIES LANGUAGE HIP)
   target_compile_options(dwarf-cloudsc-scc-parametrise-gpu-lib PRIVATE --offload-arch=gfx90a) 


    # cloudsc_add_compile_options(
    #      SOURCES dwarf_cloudsc.F90
    #      FLAGS   "-cuda") # -gpu=maxrregcount ...
    # target_include_directories(
    #     dwarf-cloudsc-scc-parametrise-gpu-lib
    #     PUBLIC
    #     ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    # )

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-loki-scc-parametrise-hip
        SOURCES dwarf_cloudsc.F90 loki-hip-scc-parametrise/cuf_cloudsc_driver_loki_mod.hip_parametrise.F90
        LIBS dwarf-cloudsc-scc-parametrise-gpu-lib
    )

    target_compile_definitions(dwarf-cloudsc-loki-scc-parametrise-hip PUBLIC CUF_DRIVER=1)

    set_target_properties(dwarf-cloudsc-loki-scc-parametrise-hip PROPERTIES LINKER_LANGUAGE Fortran)

    # set_target_properties( dwarf-cloudsc-loki-scc-parametrise-c PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    # target_link_libraries(dwarf-cloudsc-loki-scc-parametrise-c dwarf-cloudsc-scc-parametrise-gpu-lib cudadevrt cudart)

    # HOIST
    cloudsc_xmod( loki-hip-scc-hoist )

    loki_transform_convert(
        MODE hip-hoist FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
        # CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS
                ${COMMON_MODULE}/parkind1.F90
                ${COMMON_MODULE}/yomphyder.F90
                ${COMMON_MODULE}/yomcst.F90
                ${COMMON_MODULE}/yoethf.F90
                ${COMMON_MODULE}/yoecldp.F90
                ${COMMON_MODULE}/fcttre_mod.F90
                ${COMMON_MODULE}/fccld_mod.F90
        INCLUDES ${COMMON_INCLUDE}
        # DEFINITIONS CLOUDSC_GPU_TIMING ${CLOUDSC_DEFINE_STMT_FUNC}
        DEFINITIONS CLOUDSC_GPU_TIMING
        DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-hip-scc-hoist
        OUTPUT
            loki-hip-scc-hoist/cuf_cloudsc_driver_loki_mod.hip_hoist.F90
            loki-hip-scc-hoist/cloudsc_cuf_fc.F90 loki-hip-scc-hoist/cloudsc_cuf_c.c
            loki-hip-scc-hoist/yoethf_fc.F90 loki-hip-scc-hoist/yomcst_fc.F90
            loki-hip-scc-hoist/yoecldp_fc.F90
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )


      ecbuild_add_library(
        TARGET dwarf-cloudsc-scc-hoist-gpu-lib
        INSTALL_HEADERS LISTED
        SOURCES
            # yoecldp_c.h
            loki-hip-scc-hoist/cloudsc_cuf_c.c
            loki-hip-scc-hoist/cloudsc_cuf_fc.F90
            PUBLIC_INCLUDES
                $<INSTALL_INTERFACE:include>
                $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
            PUBLIC_LIBS
                cloudsc-common-lib
    )

   target_include_directories(dwarf-cloudsc-scc-hoist-gpu-lib PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

   set_source_files_properties(loki-hip-scc-hoist/cloudsc_cuf_c.c PROPERTIES LANGUAGE HIP)
   target_compile_options(dwarf-cloudsc-scc-hoist-gpu-lib PRIVATE --offload-arch=gfx90a)


    # cloudsc_add_compile_options(
    #      SOURCES dwarf_cloudsc.F90
    #      FLAGS   "-cuda") # -gpu=maxrregcount ...
    # target_include_directories(
    #     dwarf-cloudsc-scc-parametrise-gpu-lib
    #     PUBLIC
    #     ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    # )

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-loki-scc-hoist-hip
        SOURCES dwarf_cloudsc.F90 loki-hip-scc-hoist/cuf_cloudsc_driver_loki_mod.hip_hoist.F90
        LIBS dwarf-cloudsc-scc-hoist-gpu-lib
    )

    target_compile_definitions(dwarf-cloudsc-loki-scc-hoist-hip PUBLIC CUF_DRIVER=1)

    set_target_properties(dwarf-cloudsc-loki-scc-hoist-hip PROPERTIES LINKER_LANGUAGE Fortran)

    # set_target_properties( dwarf-cloudsc-loki-scc-parametrise-c PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    # target_link_libraries(dwarf-cloudsc-loki-scc-parametrise-c dwarf-cloudsc-scc-parametrise-gpu-lib cudadevrt cudart)
endif()

if (FALSE) 
    ######################################################
    ##  "Single Column Coalesced" (SCC) mode with stack ##
    ##   * Removes horizontal vector loops              ##
    ##   * Invokes compute kernel as `!$acc vector`     ##
    ##   * Allocates temporaries using pool allocator   ##
    ######################################################

    cloudsc_xmod( loki-scc-stack )

    loki_transform_convert(
        MODE scc-stack FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS ${COMMON_MODULE}/yomphyder.F90 ${COMMON_MODULE}/yoecldp.F90
        INCLUDES ${COMMON_INCLUDE}
        DEFINITIONS CLOUDSC_GPU_TIMING ${CLOUDSC_DEFINE_STMT_FUNC}
        DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-stack
        OUTPUT
            loki-scc-stack/cloudsc.scc_stack.F90
            loki-scc-stack/cloudsc_driver_loki_mod.scc_stack.F90
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-stack
        SOURCES
            dwarf_cloudsc.F90
            stack_mod.F90
            loki-scc-stack/cloudsc_driver_loki_mod.scc_stack.F90
            loki-scc-stack/cloudsc.scc_stack.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
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

    cloudsc_xmod( loki-scc-hoist )

    loki_transform_convert(
        MODE scc-hoist FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS ${COMMON_MODULE}/yomphyder.F90 ${COMMON_MODULE}/yoecldp.F90
        INCLUDES ${COMMON_INCLUDE}
        DEFINITIONS CLOUDSC_GPU_TIMING ${CLOUDSC_DEFINE_STMT_FUNC}
        DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-hoist
        OUTPUT
            loki-scc-hoist/cloudsc.scc_hoist.F90
            loki-scc-hoist/cloudsc_driver_loki_mod.scc_hoist.F90
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
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

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-hoist-serial
        COMMAND bin/dwarf-cloudsc-loki-scc-hoist
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )


    ####################################################
    ##  SCC CUF (CUDA Fortran)                        ##
    ##   * SCC with CUDA Fortran (CUF)                ##
    ##   * allows for 3 types/kind of transformation  ##
    ##    * parametrised (dimensions) local arrays    ##
    ##    * host side hoisted local arrays            ##
    ##    * dynamic memory allocation for local       ##
    ##      arrays (not recommended)                  ##
    ####################################################

endif()

if( HAVE_CUDA )
if (FALSE)
    # scc-cuf-parametrise
    cloudsc_xmod( loki-scc-cuf-parametrise )

    loki_transform_convert(
        MODE cuf-parametrise FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS
            ${COMMON_MODULE}/yomcst.F90
            ${COMMON_MODULE}/yomphyder.F90
            ${COMMON_MODULE}/yoethf.F90
            ${COMMON_MODULE}/yoecldp.F90
        INCLUDES ${COMMON_INCLUDE}
	# DEFINITIONS CLOUDSC_GPU_TIMING CUF_IMPORTS
	DEFINITIONS CLOUDSC_GPU_TIMING CUF_IMPORTS
	DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-cuf-parametrise
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

    target_compile_definitions(dwarf-cloudsc-loki-scc-cuf-parametrise PUBLIC CUF_DRIVER=1)
    target_link_options(dwarf-cloudsc-loki-scc-cuf-parametrise PUBLIC "-cuda")

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-cuf-parametrise-serial
        COMMAND bin/dwarf-cloudsc-loki-scc-cuf-parametrise
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )

    # scc-cuf-hoist
    cloudsc_xmod( loki-scc-cuf-hoist )

    loki_transform_convert(
        MODE cuf-hoist FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS
            ${COMMON_MODULE}/yomcst.F90
            ${COMMON_MODULE}/yomphyder.F90
            ${COMMON_MODULE}/yoethf.F90
            ${COMMON_MODULE}/yoecldp.F90
        INCLUDES ${COMMON_INCLUDE}
	# DEFINITIONS CLOUDSC_GPU_TIMING CUF_IMPORTS
	DEFINITIONS CLOUDSC_GPU_TIMING CUF_IMPORTS
	DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-cuf-hoist
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

    target_compile_definitions(dwarf-cloudsc-loki-scc-cuf-hoist PUBLIC CUF_DRIVER=1)
    target_link_options(dwarf-cloudsc-loki-scc-cuf-hoist PUBLIC "-cuda")

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-cuf-hoist-serial
        COMMAND bin/dwarf-cloudsc-loki-scc-cuf-hoist
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
endif()
endif()

    #############################################################################
    ## C-transpilation mode for generating vectorized C host code (experimental!)
    #############################################################################

    if ( TRUE ) # NOT HAVE_SINGLE_PRECISION )

        cloudsc_xmod( loki-c )
	
	loki_transform(
	    COMMAND convert
	    MODE c FRONTEND ${LOKI_FRONTEND} CPP
	    # CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
	    CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
	    SOURCES ${CMAKE_CURRENT_SOURCE_DIR}
	    HEADERS
	        ${COMMON_MODULE}/parkind1.F90
		${COMMON_MODULE}/yomphyder.F90
		${COMMON_MODULE}/yomcst.F90
		${COMMON_MODULE}/yoethf.F90
		${COMMON_MODULE}/yoecldp.F90
		${COMMON_MODULE}/fcttre_mod.F90
		${COMMON_MODULE}/fccld_mod.F90
	    INCLUDES ${COMMON_INCLUDE}
	    XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
	    BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-c
	    OUTPUT
	    	# loki-c/cuf_cloudsc_driver_loki_mod.c.F90
		# loki-c/cloudsc_cuf_fc.F90 loki-c/cloudsc_cuf_c.c
                loki-c/cloudsc_driver_loki_mod.c.F90
                loki-c/cloudsc_fc.F90 loki-c/cloudsc_c.c
		loki-c/yoethf_fc.F90 loki-c/yomcst_fc.F90
		loki-c/yoecldp_fc.F90
	    # DEPENDS cuf_cloudsc.F90 cuf_cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
	    DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
        )

        # Define the CLAW-CPU build target for this variant
        ecbuild_add_executable( TARGET dwarf-cloudsc-loki-c
            SOURCES
                dwarf_cloudsc.F90
		# loki-c/cuf_cloudsc_driver_loki_mod.c.F90
		# loki-c/cloudsc_cuf_fc.F90
		# loki-c/cloudsc_cuf_c.c
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
