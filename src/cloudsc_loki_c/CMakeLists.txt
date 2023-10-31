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


    #############################################################################
    ## C-transpilation mode for generating vectorized C host code (experimental!)
    #############################################################################

    # if ( NOT HAVE_SINGLE_PRECISION )
    if ( FALSE ) 
        cloudsc_xmod( loki-c )

        loki_transform_transpile(
            FRONTEND ${LOKI_FRONTEND} CPP
            HEADERS
                ${COMMON_MODULE}/parkind1.F90
                ${COMMON_MODULE}/yomphyder.F90
                ${COMMON_MODULE}/yomcst.F90
                ${COMMON_MODULE}/yoethf.F90
                ${COMMON_MODULE}/yoecldp.F90
                ${COMMON_MODULE}/fcttre_mod.F90
                ${COMMON_MODULE}/fccld_mod.F90
            DRIVER ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_driver_loki_mod.F90
	    SOURCE ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc.F90
            INCLUDES ${COMMON_INCLUDE}
            XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
            OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-c
            OUTPUT
                loki-c/cloudsc_driver_loki_mod.c.F90
                loki-c/cloudsc_fc.F90 loki-c/cloudsc_c.c
                loki-c/yoethf_fc.F90 loki-c/yomcst_fc.F90
                loki-c/yoecldp_fc.F90
            DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
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
            CONDITION HAVE_MPI
        )

    endif()

    ##
    if ( TRUE  )
    cloudsc_xmod( loki-c-scc )

    loki_transform_convert(
        MODE c FRONTEND ${LOKI_FRONTEND} CPP
	CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
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
	OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-c-scc
        OUTPUT
            loki-c-scc/cloudsc_driver_loki_mod.c.F90
            loki-c-scc/cloudsc_fc.F90 loki-c-scc/cloudsc_c.c
            loki-c-scc/yoethf_fc.F90 loki-c-scc/yomcst_fc.F90
            loki-c-scc/yoecldp_fc.F90
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )


#     loki_transform(
#	COMMAND convert
#        MODE c FRONTEND ${LOKI_FRONTEND} CPP
#        CONFIG ${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config
#	SOURCES ${CMAKE_CURRENT_SOURCE_DIR} ${COMMON_MODULE}
#	# HEADERS
#	# 	${COMMON_MODULE}/parkind1.F90
#	#         ${COMMON_MODULE}/yomphyder.F90
#	#         ${COMMON_MODULE}/yomcst.F90
#	#         ${COMMON_MODULE}/yoethf.F90
#	#         ${COMMON_MODULE}/yoecldp.F90
#	#         ${COMMON_MODULE}/fcttre_mod.F90
#	#         ${COMMON_MODULE}/fccld_mod.F90
#	INCLUDES ${COMMON_INCLUDE}
#	# DEFINITIONS CLOUDSC_GPU_TIMING ${CLOUDSC_DEFINE_STMT_FUNC}
#        DEFINITIONS CLOUDSC_GPU_TIMING
#	DATA_OFFLOAD REMOVE_OPENMP
#        XMOD ${_TARGET_XMOD_DIR} ${XMOD_DIR}
#	BUILDDIR ${CMAKE_CURRENT_BINARY_DIR}/loki-c-scc
#	OUTPUT
#            loki-c-scc/cloudsc_driver_loki_mod.F90
#            loki-c-scc/cloudsc_fc.F90 loki-c-scc/cloudsc_c.c
#            loki-c-scc/yoethf_fc.F90 loki-c-scc/yomcst_fc.F90
#            loki-c-scc/yoecldp_fc.F90
#	DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
#    )


    # ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc
    #     SOURCES
    #         dwarf_cloudsc.F90
    #         loki-scc/cloudsc_driver_loki_mod.scc.F90
    #         loki-scc/cloudsc.scc.F90
    #     LIBS
    #         cloudsc-common-lib
    #     DEFINITIONS ${CLOUDSC_DEFINITIONS}
    # )

    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-scc-c
        SOURCES
	    dwarf_cloudsc.F90
            loki-c-scc/cloudsc_driver_loki_mod.c.F90
            loki-c-scc/cloudsc_fc.F90
            loki-c-scc/cloudsc_c.c
            loki-c-scc/yoethf_fc.F90
            loki-c-scc/yomcst_fc.F90
            loki-c-scc/yoecldp_fc.F90
	    # dwarf_cloudsc.F90
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
        ENVIRONMENT "NVCOMPILER_ACC_CUDA_HEAPSIZE=64M"
    )
    endif()
    ##

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
