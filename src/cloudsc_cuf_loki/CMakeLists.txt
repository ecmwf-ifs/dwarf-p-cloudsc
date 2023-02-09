# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature

set( COMMON_MODULE "${CMAKE_CURRENT_SOURCE_DIR}/../common/module" )

ecbuild_add_option( FEATURE CLOUDSC_CUF_LOKI
    DESCRIPTION "Use Loki source-to-source transformations with CLOUDSC " DEFAULT OFF
    CONDITION HAVE_CUDA # AND HAVE_CLOUDSC_LOKI AND (Serialbox_FOUND OR HDF5_FOUND)
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

    ####################################################
    ##  SCC CUF (CUDA Fortran)                        ##
    ##   * SCC with CUDA Fortran (CUF)                ##
    ##   * allows for 3 types/kind of transformation  ##
    ##    * parametrised (dimensions) local arrays    ##
    ##    * host side hoisted local arrays            ##
    ##    * dynamic memory allocation for local       ##
    ##      arrays (not recommended)                  ##
    ####################################################

if ( HAVE_CLOUDSC_CUF_LOKI )

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
        DEFINITIONS CLOUDSC_GPU_TIMING
        DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-cuf-parametrise
        OUTPUT
            loki-scc-cuf-parametrise/cuf_cloudsc_driver_loki_mod.cuf_parametrise.F90
            loki-scc-cuf-parametrise/cuf_cloudsc.cuf_parametrise.F90
	# DEPENDS ${_OMNI_DEPENDENCIES}
	DEPENDS cuf_cloudsc.F90 cuf_cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )

    set_source_files_properties(
    loki-scc-cuf-parametrise/cuf_cloudsc_driver_loki_mod.cuf_parametrise.F90 loki-scc-cuf-parametrise/cuf_cloudsc.cuf_parametrise.F90 PROPERTIES COMPILE_FLAGS "-Mcuda=maxregcount:128")

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-cuf-parametrise
        SOURCES
            dwarf_cloudsc.F90
            loki-scc-cuf-parametrise/cuf_cloudsc_driver_loki_mod.cuf_parametrise.F90
            loki-scc-cuf-parametrise/cuf_cloudsc.cuf_parametrise.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_CUF
    )

    # target_compile_definitions(dwarf-cloudsc-loki-scc-cuf-parametrise PUBLIC USE_CUDA_DRIVER=1)
    target_link_options(dwarf-cloudsc-loki-scc-cuf-parametrise PUBLIC "-Mcuda")

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
        DEFINITIONS CLOUDSC_GPU_TIMING
        DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-cuf-hoist
        OUTPUT
            loki-scc-cuf-hoist/cuf_cloudsc_driver_loki_mod.cuf_hoist.F90
            loki-scc-cuf-hoist/cuf_cloudsc.cuf_hoist.F90
	# DEPENDS ${_OMNI_DEPENDENCIES}
	DEPENDS cuf_cloudsc.F90 cuf_cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )

    set_source_files_properties(
    loki-scc-cuf-hoist/cuf_cloudsc_driver_loki_mod.cuf_hoist.F90 loki-scc-cuf-hoist/cuf_cloudsc.cuf_hoist.F90 PROPERTIES COMPILE_FLAGS "-Mcuda=maxregcount:128")

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-cuf-hoist
        SOURCES
            dwarf_cloudsc.F90
            loki-scc-cuf-hoist/cuf_cloudsc_driver_loki_mod.cuf_hoist.F90
            loki-scc-cuf-hoist/cuf_cloudsc.cuf_hoist.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_CUF
    )

    # target_compile_definitions(dwarf-cloudsc-loki-scc-cuf-hoist PUBLIC USE_CUDA_DRIVER=1)
    target_link_options(dwarf-cloudsc-loki-scc-cuf-hoist PUBLIC "-Mcuda")

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-cuf-hoist-serial
        COMMAND bin/dwarf-cloudsc-loki-scc-cuf-hoist
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )

    # scc-cuf-dynamic
    cloudsc_xmod( loki-scc-cuf-dynamic )

    loki_transform_convert(
        MODE cuf-dynamic FRONTEND ${LOKI_FRONTEND} CPP
        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_cuf_loki.config
        PATH ${CMAKE_CURRENT_SOURCE_DIR}
        HEADERS
            ${COMMON_MODULE}/yomcst.F90
            ${COMMON_MODULE}/yomphyder.F90
            ${COMMON_MODULE}/yoethf.F90
            ${COMMON_MODULE}/yoecldp.F90
        INCLUDES ${COMMON_INCLUDE}
        DEFINITIONS CLOUDSC_GPU_TIMING
        DATA_OFFLOAD REMOVE_OPENMP
        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
        OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-scc-cuf-dynamic
        OUTPUT
            loki-scc-cuf-dynamic/cuf_cloudsc_driver_loki_mod.cuf_dynamic.F90
            loki-scc-cuf-dynamic/cuf_cloudsc.cuf_dynamic.F90
	# DEPENDS ${_OMNI_DEPENDENCIES}
	DEPENDS cuf_cloudsc.F90 cuf_cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )

    set_source_files_properties(
    loki-scc-cuf-dynamic/cuf_cloudsc_driver_loki_mod.cuf_dynamic.F90 loki-scc-cuf-dynamic/cuf_cloudsc.cuf_dynamic.F90 PROPERTIES COMPILE_FLAGS "-Mcuda=maxregcount:128")

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-cuf-dynamic
        SOURCES
            dwarf_cloudsc.F90
            loki-scc-cuf-dynamic/cuf_cloudsc_driver_loki_mod.cuf_dynamic.F90
            loki-scc-cuf-dynamic/cuf_cloudsc.cuf_dynamic.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_CUF
    )

    # target_compile_definitions(dwarf-cloudsc-loki-scc-cuf-dynamic PUBLIC USE_CUDA_DRIVER=1)
    target_link_options(dwarf-cloudsc-loki-scc-cuf-dynamic PUBLIC "-Mcuda")

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-scc-cuf-dynamic-serial
        COMMAND bin/dwarf-cloudsc-loki-scc-cuf-dynamic
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )

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
