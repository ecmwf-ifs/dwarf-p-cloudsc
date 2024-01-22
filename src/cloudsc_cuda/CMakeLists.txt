# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_C_CUDA
    DESCRIPTION "Build the CUDA version of CLOUDSC C using Serialbox" DEFAULT ON
    CONDITION Serialbox_FOUND AND HAVE_CUDA
)

if( HAVE_CLOUDSC_C_CUDA )

    enable_language(CUDA)
    enable_language(CXX)

    ###### SCC-CUDA ####
    ecbuild_add_library(
        TARGET dwarf-cloudsc-c-cuda-lib
        INSTALL_HEADERS LISTED
        SOURCES
        cloudsc/yoecldp_c.h
        cloudsc/load_state.h
            cloudsc/load_state.cu
            cloudsc/cloudsc_c.h
            cloudsc/cloudsc_c.cu
            cloudsc/cloudsc_driver.h
            cloudsc/cloudsc_driver.cu
            cloudsc/cloudsc_validate.h
            cloudsc/cloudsc_validate.cu
            cloudsc/mycpu.h
            cloudsc/mycpu.cu
        PUBLIC_INCLUDES
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>
        PUBLIC_LIBS
            Serialbox::Serialbox_C
            $<${HAVE_OMP}:OpenMP::OpenMP_C>
    )

    target_include_directories(
        dwarf-cloudsc-c-cuda-lib
        PUBLIC
        ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    )
    if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        target_compile_options(dwarf-cloudsc-c-cuda-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>>)
    else()
        target_compile_options(dwarf-cloudsc-c-cuda-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
            -gencode arch=compute_${CMAKE_CUDA_ARCHITECTURES},code=sm_${CMAKE_CUDA_ARCHITECTURES}>)
    endif()

    set_target_properties( dwarf-cloudsc-c-cuda-lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-c-cuda
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-c-cuda-lib
    )

    target_link_libraries(dwarf-cloudsc-c-cuda dwarf-cloudsc-c-cuda-lib)

    ecbuild_add_test(
        TARGET dwarf-cloudsc-c-cuda-serial
        COMMAND bin/dwarf-cloudsc-c-cuda
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ###

    ###### SCC-CUDA-HOIST ####
    ecbuild_add_library(
        TARGET dwarf-cloudsc-c-cuda-hoist-lib
        INSTALL_HEADERS LISTED
        SOURCES
            cloudsc/yoecldp_c.h
            cloudsc/load_state.h
            cloudsc/load_state.cu
            cloudsc/cloudsc_c_hoist.h
            cloudsc/cloudsc_c_hoist.cu
            cloudsc/cloudsc_driver_hoist.h
            cloudsc/cloudsc_driver_hoist.cu
            cloudsc/cloudsc_validate.h
            cloudsc/cloudsc_validate.cu
            cloudsc/mycpu.h
            cloudsc/mycpu.cu
        PUBLIC_INCLUDES
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>
        PUBLIC_LIBS
            Serialbox::Serialbox_C
            $<${HAVE_OMP}:OpenMP::OpenMP_C>
    )

    target_include_directories(
      dwarf-cloudsc-c-cuda-hoist-lib
      PUBLIC
      ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    )
    if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
      target_compile_options(dwarf-cloudsc-c-cuda-hoist-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>>)
    else()
      target_compile_options(dwarf-cloudsc-c-cuda-hoist-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
                       -gencode arch=compute_${CMAKE_CUDA_ARCHITECTURES},code=sm_${CMAKE_CUDA_ARCHITECTURES}>)
    endif()

    set_target_properties( dwarf-cloudsc-c-cuda-hoist-lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-c-cuda-hoist
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-c-cuda-hoist-lib
    )

    target_link_libraries(dwarf-cloudsc-c-cuda-hoist dwarf-cloudsc-c-cuda-hoist-lib)

    ecbuild_add_test(
        TARGET dwarf-cloudsc-c-cuda-hoist-serial
        COMMAND bin/dwarf-cloudsc-c-cuda-hoist
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ###

    ###### SCC-CUDA-K-CACHING ####
    ecbuild_add_library(
        TARGET dwarf-cloudsc-c-cuda-k-caching-lib
        INSTALL_HEADERS LISTED
        SOURCES
            cloudsc/yoecldp_c.h
            cloudsc/load_state.h
            cloudsc/load_state.cu
            cloudsc/cloudsc_c_k_caching.h
            cloudsc/cloudsc_c_k_caching.cu
            cloudsc/cloudsc_driver.h
            cloudsc/cloudsc_driver.cu
            cloudsc/cloudsc_validate.h
            cloudsc/cloudsc_validate.cu
            cloudsc/mycpu.h
            cloudsc/mycpu.cu
        PUBLIC_INCLUDES
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>
        PUBLIC_LIBS
            Serialbox::Serialbox_C
            $<${HAVE_OMP}:OpenMP::OpenMP_C>
    )

    target_include_directories(
      dwarf-cloudsc-c-cuda-k-caching-lib
      PUBLIC
      ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    )
    if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
      target_compile_options(dwarf-cloudsc-c-cuda-k-caching-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>>)
    else()
      target_compile_options(dwarf-cloudsc-c-cuda-k-caching-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
                       -gencode arch=compute_${CMAKE_CUDA_ARCHITECTURES},code=sm_${CMAKE_CUDA_ARCHITECTURES}>)
    endif()
    set_target_properties( dwarf-cloudsc-c-cuda-k-caching-lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-c-cuda-k-caching
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-c-cuda-k-caching-lib
    )

    target_link_libraries(dwarf-cloudsc-c-cuda-k-caching dwarf-cloudsc-c-cuda-k-caching-lib)

    ecbuild_add_test(
        TARGET dwarf-cloudsc-c-cuda-k-caching-serial
        COMMAND bin/dwarf-cloudsc-c-cuda-k-caching
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ###

    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )

else()
    ecbuild_info( "Serialbox and/or CUDA not found, disabling CUDA prototype(s)" )
endif()
