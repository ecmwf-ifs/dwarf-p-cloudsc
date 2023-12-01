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

    enable_language(CUDA)
    set(CMAKE_CUDA_FLAGS "-O3 -ffast-math") #  -rdc=true") #  -dc")
    find_package(CUDAToolkit)
    link_directories("/usr/local/apps/nvidia/22.11/Linux_x86_64/22.11/math_libs/lib64")

    if (TRUE) 
    #### CUDA PARAMETRISE #####
    cloudsc_xmod( loki-c-scc-parametrise )

    loki_transform_convert(
        MODE c-parametrise FRONTEND ${LOKI_FRONTEND} CPP
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
	OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-c-scc-parametrise
        OUTPUT
            loki-c-scc-parametrise/cloudsc_driver_loki_mod.c_parametrise.F90
            loki-c-scc-parametrise/cloudsc_fc.F90 loki-c-scc-parametrise/cloudsc_c.c
            loki-c-scc-parametrise/yoethf_fc.F90 loki-c-scc-parametrise/yomcst_fc.F90
            loki-c-scc-parametrise/yoecldp_fc.F90 loki-c-scc-parametrise/cloudsc_init_seq_c.c 
	    loki-c-scc-parametrise/cloudsc_init_vec_c.c
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )


      ecbuild_add_library(
        TARGET dwarf-cloudsc-scc-parametrise-modified-gpu-lib
        INSTALL_HEADERS LISTED
        SOURCES
	    yoecldp_c.h
            loki-c-scc-parametrise/cloudsc_c.c
    	    loki-c-scc-parametrise/cloudsc_fc.F90
	    # loki-c-scc-parametrise/cloudsc_init_seq_c.c
	    # loki-c-scc-parametrise/cloudsc_init_vec_c.c
	    PUBLIC_INCLUDES
    	    	$<INSTALL_INTERFACE:include>
            	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
    	    PUBLIC_LIBS
    	        cloudsc-common-lib
    )

    target_include_directories(dwarf-cloudsc-scc-parametrise-modified-gpu-lib PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

    set_source_files_properties(loki-c-scc-parametrise/cloudsc_c.c loki-c-scc-parametrise/cloudsc_init_seq_c.c loki-c-scc-parametrise/cloudsc_init_vec_c.c PROPERTIES LANGUAGE CUDA)

   target_compile_options(dwarf-cloudsc-scc-parametrise-modified-gpu-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
	   --generate-line-info
	   >) 


      cloudsc_add_compile_options(
	 SOURCES dwarf_cloudsc.F90
         FLAGS   "-cuda") # -gpu=maxrregcount ...
    target_include_directories(
        dwarf-cloudsc-scc-parametrise-modified-gpu-lib
        PUBLIC
        ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    )

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-loki-scc-parametrise-modified-c
	SOURCES dwarf_cloudsc.F90 loki-c-scc-parametrise/cloudsc_driver_loki_mod.c_parametrise.F90
	LIBS dwarf-cloudsc-scc-parametrise-modified-gpu-lib
    )

    set_target_properties(dwarf-cloudsc-loki-scc-parametrise-modified-c PROPERTIES LINKER_LANGUAGE Fortran)
    # target_link_options(dwarf-cloudsc-loki-scc-parametrise-modified-c PUBLIC "-cuda")
    # set_target_properties( dwarf-cloudsc-loki-scc-parametrise-modified-c PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_link_libraries(dwarf-cloudsc-loki-scc-parametrise-modified-c dwarf-cloudsc-scc-parametrise-modified-gpu-lib cudadevrt cudart)
    #### 
    endif()


    #### CUDA HOIST #####
    cloudsc_xmod( loki-c-scc-hoist )

    loki_transform_convert(
        MODE c-hoist FRONTEND ${LOKI_FRONTEND} CPP
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
	OUTPATH ${CMAKE_CURRENT_BINARY_DIR}/loki-c-scc-hoist
        OUTPUT
            loki-c-scc-hoist/cloudsc_driver_loki_mod.c_hoist.F90
            loki-c-scc-hoist/cloudsc_fc.F90 loki-c-scc-hoist/cloudsc_c.c
            loki-c-scc-hoist/yoethf_fc.F90 loki-c-scc-hoist/yomcst_fc.F90
            loki-c-scc-hoist/yoecldp_fc.F90 loki-c-scc-hoist/cloudsc_init_seq_c.c 
	    loki-c-scc-hoist/cloudsc_init_vec_c.c
        DEPENDS cloudsc.F90 cloudsc_driver_loki_mod.F90 ${_OMNI_DEPENDENCIES}
    )

     #  ecbuild_add_library(
     #	TARGET dwarf-cloudsc-scc-hoist-gpu-device-lib
     #	SOURCES
     #       loki-c-scc-hoist/cloudsc_init_seq_c.c
     #       loki-c-scc-hoist/cloudsc_init_vec_c.c
     #       PUBLIC_INCLUDES
     #           $<INSTALL_INTERFACE:include>
     #           $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
     #	TYPE STATIC
     # )
      
     # set_source_files_properties(loki-c-scc-hoist/cloudsc_init_seq_c.c loki-c-scc-hoist/cloudsc_init_vec_c.c PROPERTIES LANGUAGE CUDA)

      ecbuild_add_library(
        TARGET dwarf-cloudsc-scc-hoist-modified-gpu-lib
        INSTALL_HEADERS LISTED
        SOURCES
	    yoecldp_c.h
            loki-c-scc-hoist/cloudsc_c.c
	    # loki-c-scc-hoist/cloudsc_init_seq_c.c
	    # loki-c-scc-hoist/cloudsc_init_vec_c.c
    	    loki-c-scc-hoist/cloudsc_fc.F90
	    PUBLIC_INCLUDES
    	    	$<INSTALL_INTERFACE:include>
            	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
    	    PUBLIC_LIBS
    	        cloudsc-common-lib
		# TYPE STATIC
    )

    # set_target_properties( dwarf-cloudsc-scc-hoist-gpu-lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    # target_link_libraries(dwarf-cloudsc-scc-hoist-gpu-lib PUBLIC dwarf-cloudsc-scc-hoist-gpu-device-lib)    

    target_include_directories(dwarf-cloudsc-scc-hoist-modified-gpu-lib PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

    set_source_files_properties(loki-c-scc-hoist/cloudsc_c.c loki-c-scc-hoist/cloudsc_init_seq_c.c loki-c-scc-hoist/cloudsc_init_vec_c.c PROPERTIES LANGUAGE CUDA)
    # set_source_files_properties(loki-c-scc-hoist/cloudsc_c.c PROPERTIES LANGUAGE CUDA)

    target_compile_options(dwarf-cloudsc-scc-hoist-modified-gpu-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
	   --generate-line-info
	   >) 


      cloudsc_add_compile_options(
	 SOURCES dwarf_cloudsc.F90
         FLAGS   "-cuda") # -gpu=maxrregcount ...
    target_include_directories(
        dwarf-cloudsc-scc-hoist-modified-gpu-lib
        PUBLIC
        ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    )

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-loki-scc-hoist-modified-c
	SOURCES dwarf_cloudsc.F90 loki-c-scc-hoist/cloudsc_driver_loki_mod.c_hoist.F90
	LIBS dwarf-cloudsc-scc-hoist-modified-gpu-lib
    )

    set_target_properties(dwarf-cloudsc-loki-scc-hoist-modified-c PROPERTIES LINKER_LANGUAGE Fortran)
    # target_link_options(dwarf-cloudsc-loki-scc-hoist-modified-c PUBLIC "-cuda")    
    ## set_target_properties( dwarf-cloudsc-loki-scc-hoist-c PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    target_link_libraries(dwarf-cloudsc-loki-scc-hoist-modified-c dwarf-cloudsc-scc-hoist-modified-gpu-lib cudadevrt cudart)
    #### 

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
