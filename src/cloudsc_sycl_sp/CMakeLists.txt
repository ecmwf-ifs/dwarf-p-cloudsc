# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_SYCL
    DESCRIPTION "Build the SYCL version CLOUDSC using Serialbox" DEFAULT ON
    CONDITION (Serialbox_FOUND OR HDF5_FOUND) AND HAVE_SYCL
)

if( HAVE_CLOUDSC_SYCL AND HAVE_SINGLE_PRECISION )

    unset(SYCL_COMPILE_OPTIONS)
    if( CMAKE_Fortran_COMPILER_ID MATCHES "IntelLLVM" )
    set(SYCL_COMPILE_OPTIONS -gline-tables-only -fsycl -fsycl-targets=spir64 -fsycl-early-optimizations -Xsycl-target-backend "-options -cl-fast-relaxed-match")
    endif()

    if (CLOUDSC_MATH STREQUAL "C")
        set(GPU_MATH "1")
        ecbuild_info("Building SYCL variant using C math interface")
    elseif (CLOUDSC_MATH STREQUAL "STD")
        set(GPU_MATH "2")
        ecbuild_info("Building SYCL variant using C++ STL math interface")
    elseif(CLOUDSC_MATH STREQUAL "SYCL")
        set(GPU_MATH "3")
        ecbuild_info("Building SYCL variant using SYCL math interface")
    else ()
        set(GPU_MATH "1")
        ecbuild_info("Building SYCL variant using C math interface")
    endif()

    if (SYCL_SUB_GROUP_SIZE)
        set(SUB_GROUP_SIZE "${SYCL_SUB_GROUP_SIZE}")
    else ()
        set(SUB_GROUP_SIZE "16")
    endif()

    ecbuild_add_library(
        TARGET dwarf-cloudsc-scc-sycl-lib
        INSTALL_HEADERS LISTED
        SOURCES
            cloudsc/yoecldp_c.h
            cloudsc/load_state.h
            cloudsc/load_state.cpp
            cloudsc/cloudsc_c.kernel
            cloudsc/cloudsc_driver.h
            cloudsc/cloudsc_driver.cpp
            cloudsc/cloudsc_validate.h
            cloudsc/cloudsc_validate.cpp
            cloudsc/mycpu.h
            cloudsc/mycpu.cpp
        PUBLIC_INCLUDES
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>
        PUBLIC_LIBS
            $<${HAVE_HDF5}:hdf5::hdf5>
            $<${HAVE_SERIALBOX}:Serialbox::Serialbox_C>
        $<${HAVE_OMP}:OpenMP::OpenMP_C>
    DEFINITIONS
        ${CLOUDSC_DEFINITIONS}
    )

    target_compile_definitions(dwarf-cloudsc-scc-sycl-lib PUBLIC C_MATH=1)
    target_compile_definitions(dwarf-cloudsc-scc-sycl-lib PUBLIC STD_MATH=2)
    target_compile_definitions(dwarf-cloudsc-scc-sycl-lib PUBLIC SYCL_MATH=3)
    target_compile_definitions(dwarf-cloudsc-scc-sycl-lib PUBLIC GPU_MATH=${GPU_MATH})
    target_compile_definitions(dwarf-cloudsc-scc-sycl-lib PUBLIC SUB_GROUP_SIZE=${SUB_GROUP_SIZE})

    if( DEFINED SYCL_COMPILE_OPTIONS )
        target_compile_options(dwarf-cloudsc-scc-sycl-lib PUBLIC ${SYCL_COMPILE_OPTIONS})
    endif()

    add_sycl_to_target(
        TARGET
            dwarf-cloudsc-scc-sycl-lib
        SOURCES
            cloudsc/cloudsc_c.kernel
    )

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-scc-sycl
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-scc-sycl-lib
    )

    if( DEFINED SYCL_COMPILE_OPTIONS )
        target_compile_options(dwarf-cloudsc-scc-sycl PUBLIC ${SYCL_COMPILE_OPTIONS})
    endif()

    ecbuild_add_test(
        TARGET dwarf-cloudsc-sycl-serial
        COMMAND bin/dwarf-cloudsc-scc-sycl
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )

    ######
    ecbuild_add_library(
        TARGET dwarf-cloudsc-scc-hoist-sycl-lib
        INSTALL_HEADERS LISTED
        SOURCES
            cloudsc/yoecldp_c.h
            cloudsc/load_state.h
            cloudsc/load_state.cpp
            cloudsc/cloudsc_c_hoist.kernel
            cloudsc/cloudsc_driver_hoist.h
            cloudsc/cloudsc_driver_hoist.cpp
            cloudsc/cloudsc_validate.h
            cloudsc/cloudsc_validate.cpp
            cloudsc/mycpu.h
            cloudsc/mycpu.cpp
        PUBLIC_INCLUDES
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>
        PUBLIC_LIBS
            $<${HAVE_HDF5}:hdf5::hdf5>
            $<${HAVE_SERIALBOX}:Serialbox::Serialbox_C>
            $<${HAVE_OMP}:OpenMP::OpenMP_C>
        DEFINITIONS
            ${CLOUDSC_DEFINITIONS}
    )

    target_compile_definitions(dwarf-cloudsc-scc-hoist-sycl-lib PUBLIC C_MATH=1)
    target_compile_definitions(dwarf-cloudsc-scc-hoist-sycl-lib PUBLIC STD_MATH=2)
    target_compile_definitions(dwarf-cloudsc-scc-hoist-sycl-lib PUBLIC SYCL_MATH=3)
    target_compile_definitions(dwarf-cloudsc-scc-hoist-sycl-lib PUBLIC GPU_MATH=${GPU_MATH})
    target_compile_definitions(dwarf-cloudsc-scc-hoist-sycl-lib PUBLIC SUB_GROUP_SIZE=${SUB_GROUP_SIZE})

    if( DEFINED SYCL_COMPILE_OPTIONS )
        target_compile_options(dwarf-cloudsc-scc-hoist-sycl-lib PUBLIC ${SYCL_COMPILE_OPTIONS})
    endif()

    add_sycl_to_target(
        TARGET
            dwarf-cloudsc-scc-hoist-sycl-lib
        SOURCES
            cloudsc/cloudsc_c_hoist.kernel
    )

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-scc-hoist-sycl
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-scc-hoist-sycl-lib
    )
    if( DEFINED SYCL_COMPILE_OPTIONS )
        target_compile_options(dwarf-cloudsc-scc-hoist-sycl PUBLIC ${SYCL_COMPILE_OPTIONS})
    endif()

    ecbuild_add_test(
        TARGET dwarf-cloudsc-sycl-hoist-serial
        COMMAND bin/dwarf-cloudsc-scc-hoist-sycl
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )

    ######

    ecbuild_add_library(
        TARGET dwarf-cloudsc-scc-k-caching-sycl-lib
        INSTALL_HEADERS LISTED
        SOURCES
            cloudsc/yoecldp_c.h
            cloudsc/load_state.h
            cloudsc/load_state.cpp
            cloudsc/cloudsc_c_k_caching.kernel
            cloudsc/cloudsc_driver.h
            cloudsc/cloudsc_driver_k_caching.cpp
            cloudsc/cloudsc_validate.h
            cloudsc/cloudsc_validate.cpp
            cloudsc/mycpu.h
            cloudsc/mycpu.cpp
        PUBLIC_INCLUDES
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>
        PUBLIC_LIBS
            $<${HAVE_HDF5}:hdf5::hdf5>
            $<${HAVE_SERIALBOX}:Serialbox::Serialbox_C>
            $<${HAVE_OMP}:OpenMP::OpenMP_C>
        DEFINITIONS
            ${CLOUDSC_DEFINITIONS}
    )

    target_compile_definitions(dwarf-cloudsc-scc-k-caching-sycl-lib PUBLIC C_MATH=1)
    target_compile_definitions(dwarf-cloudsc-scc-k-caching-sycl-lib PUBLIC STD_MATH=2)
    target_compile_definitions(dwarf-cloudsc-scc-k-caching-sycl-lib PUBLIC SYCL_MATH=3)
    target_compile_definitions(dwarf-cloudsc-scc-k-caching-sycl-lib PUBLIC GPU_MATH=${GPU_MATH})
    target_compile_definitions(dwarf-cloudsc-scc-k-caching-sycl-lib PUBLIC SUB_GROUP_SIZE=${SUB_GROUP_SIZE})

    if( DEFINED SYCL_COMPILE_OPTIONS )
        target_compile_options(dwarf-cloudsc-scc-k-caching-sycl-lib PUBLIC ${SYCL_COMPILE_OPTIONS})
    endif()

    add_sycl_to_target(
        TARGET
            dwarf-cloudsc-scc-k-caching-sycl-lib
        SOURCES
            cloudsc/cloudsc_c_k_caching.kernel
    )

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-scc-k-caching-sycl
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-scc-k-caching-sycl-lib
    )
    if( DEFINED SYCL_COMPILE_OPTIONS )
        target_compile_options(dwarf-cloudsc-scc-k-caching-sycl PUBLIC ${SYCL_COMPILE_OPTIONS})
    endif()


    ecbuild_add_test(
        TARGET dwarf-cloudsc-sycl-k-caching-serial
        COMMAND bin/dwarf-cloudsc-scc-k-caching-sycl
        ARGS 1 1000 128
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
