# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

if( NOT DEFINED SYCL_TARGETS )
    set( SYCL_TARGETS "spir64" )
endif()

if( CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM" )
    if( NOT DEFINED SYCL_COMPILE_OPTIONS )
        set( SYCL_COMPILE_OPTIONS -gline-tables-only -fsycl -fsycl-targets=${SYCL_TARGETS} -fsycl-early-optimizations )
    endif()
endif()

set(CLOUDSC_SYCL_DEFINITIONS
    C_MATH=1
    STD_MATH=2
    SYCL_MATH=3
)

if( NOT DEFINED CLOUDSC_CPP_MATH )
    set( CLOUDSC_CPP_MATH "C" )
endif()

if (CLOUDSC_CPP_MATH STREQUAL "C")
    list(APPEND CLOUDSC_SYCL_DEFINITIONS GPU_MATH=1)
    ecbuild_info("Building SYCL variant using C math interface")
elseif (CLOUDSC_CPP_MATH STREQUAL "STD")
    list(APPEND CLOUDSC_SYCL_DEFINITIONS GPU_MATH=2)
    ecbuild_info("Building SYCL variant using C++ STL math interface")
elseif(CLOUDSC_CPP_MATH STREQUAL "SYCL")
    list(APPEND CLOUDSC_SYCL_DEFINITIONS GPU_MATH=3)
    ecbuild_info("Building SYCL variant using SYCL math interface")
endif()

if (SYCL_SUB_GROUP_SIZE)
    list(APPEND CLOUDSC_SYCL_DEFINITIONS SUB_GROUP_SIZE=${SYCL_SUB_GROUP_SIZE})
else ()
    list(APPEND CLOUDSC_SYCL_DEFINITIONS SUB_GROUP_SIZE=16)
endif()

### Helper macro to create a target and associated tests

function( cloudsc_add_sycl_exe )

    set( options "" )
    set( single_value_args TARGET KERNEL )
    set( multi_value_args DRIVER )

    cmake_parse_arguments( _CL "${options}" "${single_value_args}" "${multi_value_args}" ${ARGN} )

    if( _CL_UNPARSED_ARGUMENTS )
        ecbuild_critical( "Unknown keywords given to cloudsc_add_loki_exe(): \"${_CL_UNPARSED_ARGUMENTS}\"")
    endif()

    ecbuild_add_library(
        TARGET ${_CL_TARGET}-${prec}-lib
        INSTALL_HEADERS LISTED
        SOURCES
            cloudsc/yoecldp_c.h
            cloudsc/dtype.h
            cloudsc/load_state.h
            cloudsc/load_state.cpp
            cloudsc/cloudsc_validate.h
            cloudsc/cloudsc_validate.cpp
            cloudsc/mycpu.h
            cloudsc/mycpu.cpp
            ${_CL_DRIVER}
            ${_CL_KERNEL}
        PUBLIC_INCLUDES
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>
        PUBLIC_LIBS
            $<${HAVE_HDF5}:hdf5::hdf5>
            $<${HAVE_SERIALBOX}:Serialbox::Serialbox_C>
            $<${HAVE_OMP}:OpenMP::OpenMP_C>
        DEFINITIONS
            ${CLOUDSC_DEFINITIONS}
            ${CLOUDSC_SYCL_DEFINITIONS}
            $<$<STREQUAL:"${prec}","sp">:SINGLE>
    )

    add_sycl_to_target(
        TARGET  ${_CL_TARGET}-${prec}-lib
        SOURCES ${_CL_DRIVER} ${_CL_KERNEL}
    )

    ecbuild_add_executable(
        TARGET      ${_CL_TARGET}-${prec}
        SOURCES     dwarf_cloudsc.cpp
        LIBS        ${_CL_TARGET}-${prec}-lib
    )

    if( DEFINED SYCL_COMPILE_OPTIONS )
        target_compile_options( ${_CL_TARGET}-${prec}-lib
            PUBLIC ${SYCL_COMPILE_OPTIONS} )
        target_link_options( ${_CL_TARGET}-${prec}-lib
            PUBLIC ${SYCL_COMPILE_OPTIONS} )
        target_compile_options( ${_CL_TARGET}-${prec}
            PUBLIC ${SYCL_COMPILE_OPTIONS} )
        target_link_options( ${_CL_TARGET}-${prec}
            PUBLIC ${SYCL_COMPILE_OPTIONS} )
    endif()

    ecbuild_add_test(
        TARGET      ${_CL_TARGET}-${prec}-serial
        COMMAND     ${_CL_TARGET}-${prec}
        ARGS        1 1024 128
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )

endfunction()

### Kernel variants

# SCC
cloudsc_add_sycl_exe(
    TARGET dwarf-cloudsc-sycl-scc
    DRIVER
        cloudsc/cloudsc_driver.h
        cloudsc/cloudsc_driver.cpp
    KERNEL
        cloudsc/cloudsc_c.kernel
)

# SCC-HOIST
cloudsc_add_sycl_exe(
    TARGET dwarf-cloudsc-sycl-scc-hoist
    DRIVER
        cloudsc/cloudsc_driver_hoist.h
        cloudsc/cloudsc_driver_hoist.cpp
    KERNEL
        cloudsc/cloudsc_c_hoist.kernel
)

# SCC-K-CACHING
cloudsc_add_sycl_exe(
    TARGET dwarf-cloudsc-sycl-scc-k-caching
    DRIVER
        cloudsc/cloudsc_driver.h
        cloudsc/cloudsc_driver_k_caching.cpp
    KERNEL
        cloudsc/cloudsc_c_k_caching.kernel
)
