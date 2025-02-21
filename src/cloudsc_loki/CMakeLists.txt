# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

set( LOKI_FRONTEND "fp" CACHE STRING "Frontend parser for Loki transforms" )
set( XMOD_DIR "${CMAKE_CURRENT_SOURCE_DIR}/xmod" )

### Select appropriate common lib for GPU targets

# Currently, directive-based GPU recipes in Loki generate only OpenACC code, therefore
# we link to the acc-variety of the common lib if it is available. In the future, if
# other backends are added to Loki (e.g., OpenMP target support), we can select the
# common lib corresponding to the transformation mode

if( HAVE_ACC )
    set( GPU_COMMON_LIB acc )
else()
    set( GPU_COMMON_LIB lib )
endif()

### Helper macro to create Loki transformation variants

# This defines a shared library with kernel and driver, applies loki_transform_target to that
# library and links this into an executable with corresponding test definitions

function( cloudsc_add_loki_exe )

    set( options CPP_PLAN )
    set( single_value_args TARGET MODE CONFIG COMMON_LIB DRIVER KERNEL )
    set( multi_value_args DEFINITIONS TEST_ENV )

    cmake_parse_arguments( _CL "${options}" "${single_value_args}" "${multi_value_args}" ${ARGN} )

    if( _CL_UNPARSED_ARGUMENTS )
        ecbuild_critical( "Unknown keywords given to cloudsc_add_loki_exe(): \"${_CL_UNPARSED_ARGUMENTS}\"")
    endif()

    if( NOT _CL_DRIVER )
        set( _CL_DRIVER cloudsc_driver_loki_mod.F90 )
    endif()
    if( NOT _CL_KERNEL )
        set( _CL_KERNEL cloudsc.F90 )
    endif()

    set( TRAFO_OPTIONS )
    if( _CL_CPP_PLAN )
        list( APPEND TRAFO_OPTIONS CPP_PLAN )
    endif()

    ecbuild_add_library(
        TARGET      ${_CL_TARGET}-${prec}-lib
        TYPE        STATIC
        SOURCES
            ${_CL_DRIVER}
            ${_CL_KERNEL}
        PUBLIC_DEFINITIONS
            ${_CL_DEFINITIONS}
        PUBLIC_INCLUDES
            "$<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/module-${_CL_TARGET}>"
        PUBLIC_LIBS
            cloudsc-common-${prec}-${_CL_COMMON_LIB}
    )

    set_target_properties( ${_CL_TARGET}-${prec}-lib
        PROPERTIES
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module-${_CL_TARGET}"
    )

    loki_transform_target(
        TARGET      ${_CL_TARGET}-${prec}-lib
        COMMAND     convert
        MODE        ${_CL_MODE}
        FRONTEND    ${LOKI_FRONTEND}
        CONFIG      "${CMAKE_CURRENT_SOURCE_DIR}/${_CL_CONFIG}"
        PLAN        "${CMAKE_CURRENT_BINARY_DIR}/plan-${_CL_MODE}.cmake"
        CPP
        DEFINITIONS ${_CL_DEFINITIONS}
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        # XMOD
        #     ${CMAKE_CURRENT_BINARY_DIR}
        #     ${XMOD_DIR}
        ${TRAFO_OPTIONS}
    )

    ecbuild_add_executable(
        TARGET ${_CL_TARGET}-${prec}
        SOURCES dwarf_cloudsc.F90
        LIBS
            ${_CL_TARGET}-${prec}-lib
            cloudsc-common-${prec}-${_CL_COMMON_LIB}
    )

    ecbuild_add_test(
        TARGET ${_CL_TARGET}-${prec}-serial
        COMMAND ${_CL_TARGET}-${prec}
        ARGS 1 100 32
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        OMP 1
        ENVIRONMENT "${_CL_TEST_ENV}"
    )
    ecbuild_add_test(
        TARGET ${_CL_TARGET}-${prec}-mpi
        COMMAND ${_CL_TARGET}-${prec}
        ARGS 1 100 32
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        MPI 2
        OMP 1
        CONDITION HAVE_MPI
        ENVIRONMENT "${_CL_TEST_ENV}"
    )

endfunction()


####################################################
##  Idempotence mode:                             ##
##   * Internal "do-nothing" mode for Loki debug  ##
####################################################
cloudsc_add_loki_exe(
    TARGET dwarf-cloudsc-loki-idem
    MODE idem
    CONFIG cloudsc_loki.config
    COMMON_LIB lib
    DEFINITIONS ${CLOUDSC_DEFINITIONS}
)

############################################################
##  Idempotence mode with pool-allocator for temporaries: ##
##   * Internal "do-nothing" mode for Loki debug,         ##
##     but with temporaries handling                      ##
############################################################
cloudsc_add_loki_exe(
    TARGET dwarf-cloudsc-loki-idem-stack
    MODE idem-stack
    CONFIG cloudsc_loki.config
    COMMON_LIB lib
    DEFINITIONS ${CLOUDSC_DEFINITIONS}
)

if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_compile_options(
        dwarf-cloudsc-loki-idem-stack-${prec}-lib
        PRIVATE "-fcray-pointer" )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
    target_compile_options(
        dwarf-cloudsc-loki-idem-stack-${prec}-lib
        PRIVATE "-Mcray=pointer" )
elseif( CMAKE_FORTRAN_COMPILER_ID MATCHES "Intel" )
    target_compile_options(
        dwarf-cloudsc-loki-idem-stack-${prec}-lib
        PRIVATE "-safe-cray-ptr" )
endif()

####################################################
##  "Single Column Coalesced" (SCC) mode          ##
##   * Removes horizontal vector loops            ##
##   * Invokes compute kernel as `!$acc vector`   ##
####################################################
cloudsc_add_loki_exe(
    TARGET dwarf-cloudsc-loki-scc
    MODE scc
    CONFIG cloudsc_loki.config
    COMMON_LIB ${GPU_COMMON_LIB}
    DEFINITIONS
        ${CLOUDSC_DEFINITIONS}
        CLOUDSC_GPU_TIMING
    TEST_ENV "NVCOMPILER_ACC_CUDA_HEAPSIZE=128M"
)

######################################################
##  "Single Column Coalesced" (SCC) mode with stack ##
##   * Removes horizontal vector loops              ##
##   * Invokes compute kernel as `!$acc vector`     ##
##   * Allocates temporaries using pool allocator   ##
######################################################
cloudsc_add_loki_exe(
    TARGET dwarf-cloudsc-loki-scc-stack
    MODE scc-stack
    CONFIG cloudsc_loki.config
    COMMON_LIB ${GPU_COMMON_LIB}
    DEFINITIONS
        ${CLOUDSC_DEFINITIONS}
        CLOUDSC_GPU_TIMING
)

if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
    target_compile_options(
        dwarf-cloudsc-loki-scc-stack-${prec}-lib
        PRIVATE "-fcray-pointer" )
elseif( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
    target_compile_options(
        dwarf-cloudsc-loki-scc-stack-${prec}-lib
        PRIVATE "-Mcray=pointer" )
endif()

####################################################
##  SCC-hoist mode                                ##
##   * SCC with vector loop hoisted               ##
##   * Kernel is "seq, but args are full blocks   ##
##   * Temporary arrays hoisted to driver         ##
####################################################
cloudsc_add_loki_exe(
    TARGET dwarf-cloudsc-loki-scc-hoist
    MODE scc-hoist
    CONFIG cloudsc_loki.config
    COMMON_LIB ${GPU_COMMON_LIB}
    DEFINITIONS
        ${CLOUDSC_DEFINITIONS}
        CLOUDSC_GPU_TIMING
)

if( HAVE_FIELD_API )
    ####################################################
    ##  "Single Column Coalesced Field" (SCC-FIELD)   ##
    ##   * Removes horizontal vector loops            ##
    ##   * Invokes compute kernel as `!$acc vector`   ##
    ##   * Uses Field API for device data offloading  ##
    ####################################################
    cloudsc_add_loki_exe(
        TARGET dwarf-cloudsc-loki-scc-field
        MODE scc-field
        CONFIG cloudsc_loki_field_offload.config
        COMMON_LIB field_api
        DEFINITIONS
            ${CLOUDSC_DEFINITIONS}
            CLOUDSC_GPU_TIMING
            CLOUDSC_GPU_SCC_FIELD
        DRIVER cloudsc_driver_field_loki_mod.F90
        TEST_ENV "NVCOMPILER_ACC_CUDA_HEAPSIZE=128M"
    )
endif()

if( HAVE_CUDA )
    ####################################################
    ##  SCC CUF (CUDA Fortran)                        ##
    ##   * SCC with CUDA Fortran (CUF)                ##
    ##   * allows for 3 types/kind of transformation  ##
    ##    * parametrised (dimensions) local arrays    ##
    ##    * host side hoisted local arrays            ##
    ##    * dynamic memory allocation for local       ##
    ##      arrays (not recommended)                  ##
    ####################################################)

    foreach( variant cuf-parametrise cuf-hoist )

        string( REPLACE "-" "_" variant_ ${variant} )

        cloudsc_add_loki_exe(
            TARGET dwarf-cloudsc-loki-gpu-scc-${variant}
            MODE ${variant}
            CONFIG cloudsc_cuf_loki.config
            COMMON_LIB cuda
            DEFINITIONS
                ${CLOUDSC_DEFINITIONS}
                CLOUDSC_GPU_TIMING
                CLOUDSC_GPU_SCC_CUF
            DRIVER cuf_cloudsc_driver_loki_mod.F90
            KERNEL cuf_cloudsc.F90
        )

        target_link_options( dwarf-cloudsc-loki-gpu-scc-${variant}-${prec}-lib
            PUBLIC "-cuda"
        )

        set_source_files_properties(
            cuf_cloudsc_driver_loki_mod.${variant_}.F90
            cuf_cloudsc.${variant_}.F90
            PROPERTIES COMPILE_FLAGS "-cuda -gpu=maxregcount:128"
        )

    endforeach()

endif()

if( "${prec}" STREQUAL "dp" )

    ####################################################
    ## C-transpilation mode for generating vectorized ##
    ## C host code (experimental!)                    ##
    ####################################################

    # TODO: This uses an outdated entry point to loki transformations
    # because the C-transpilation is not yet fully compatible with
    # Loki's planning functionality.

    loki_transform(
        COMMAND     convert
        MODE        c
        FRONTEND    ${LOKI_FRONTEND}
        CPP
        CONFIG      "${CMAKE_CURRENT_SOURCE_DIR}/cloudsc_loki.config"
        SOURCES
            ${CMAKE_CURRENT_SOURCE_DIR}
            ${COMMON_MODULE}
        INCLUDES
            ${COMMON_INCLUDE}
        XMOD        ${XMOD_DIR}
        BUILDDIR    ${CMAKE_CURRENT_BINARY_DIR}
        OUTPUT
            cloudsc_driver_loki_mod.c.F90
            cloudsc_fc.F90
            cloudsc_c.c
            yoethf_fc.F90
            yomcst_fc.F90
            yoecldp_fc.F90
        DEPENDS
            cloudsc.F90
            cloudsc_driver_loki_mod.F90
    )

    # Define the build target for the Loki-C variant
    ecbuild_add_executable( TARGET dwarf-cloudsc-loki-c-${prec}
        SOURCES
            dwarf_cloudsc.F90
            cloudsc_driver_loki_mod.c.F90
            cloudsc_fc.F90
            cloudsc_c.c
            yoethf_fc.F90
            yomcst_fc.F90
            yoecldp_fc.F90
        LIBS
            cloudsc-common-${prec}-lib
        DEFINITIONS
            ${CLOUDSC_DEFINITIONS}
    )

    set_target_properties( dwarf-cloudsc-loki-c-${prec}
        PROPERTIES
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module-dwarf-cloudsc-loki-c"
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-c-${prec}-serial
        COMMAND dwarf-cloudsc-loki-c-${prec}
        ARGS 1 100 32
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        OMP 1
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-loki-c-${prec}-omp
        COMMAND dwarf-cloudsc-loki-c-${prec}
        ARGS 4 100 32
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        OMP 4
        CONDITION HAVE_OMP
    )

endif()
