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
  if ( HAVE_SINGLE_PRECISION )
    enable_language(CUDA)
    enable_language(CXX)
    set(CMAKE_CUDA_FLAGS "${CMAKE_CUDA_FLAGS} --ptxas-options=-v -O0 -g")

   if (TRUE) 
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
            --ptxas-options=-v -O3 -use_fast_math -maxrregcount=128 -gencode arch=compute_${CMAKE_CUDA_ARCHITECTURES},code=sm_${CMAKE_CUDA_ARCHITECTURES}>)
    endif()

    set_target_properties( dwarf-cloudsc-c-cuda-lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-c-cuda-sp
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-c-cuda-lib
    )

    target_link_libraries(dwarf-cloudsc-c-cuda-sp dwarf-cloudsc-c-cuda-lib)

    ecbuild_add_test(
        TARGET dwarf-cloudsc-c-cuda-sp-serial
        COMMAND bin/dwarf-cloudsc-c-cuda-sp
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ###
endif()

if (TRUE)
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
                       --ptxas-options=-v -O3 -use_fast_math -maxrregcount=128 -gencode arch=compute_${CMAKE_CUDA_ARCHITECTURES},code=sm_${CMAKE_CUDA_ARCHITECTURES}>)
    endif()

    set_target_properties( dwarf-cloudsc-c-cuda-hoist-lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-c-cuda-hoist-sp
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-c-cuda-hoist-lib
    )

    target_link_libraries(dwarf-cloudsc-c-cuda-hoist-sp dwarf-cloudsc-c-cuda-hoist-lib)

    ecbuild_add_test(
        TARGET dwarf-cloudsc-c-cuda-hoist-sp-serial
        COMMAND bin/dwarf-cloudsc-c-cuda-hoist-sp
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ###
endif()
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
            cloudsc/cloudsc_driver_k_caching.h
            cloudsc/cloudsc_driver_k_caching.cu
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
                       --ptxas-options=-v -O3 -use_fast_math -maxrregcount=128 -gencode arch=compute_${CMAKE_CUDA_ARCHITECTURES},code=sm_${CMAKE_CUDA_ARCHITECTURES}>)
    endif()
    set_target_properties( dwarf-cloudsc-c-cuda-k-caching-lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-c-cuda-k-caching-sp
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-c-cuda-k-caching-lib
    )

    target_link_libraries(dwarf-cloudsc-c-cuda-k-caching-sp dwarf-cloudsc-c-cuda-k-caching-lib)

    ecbuild_add_test(
        TARGET dwarf-cloudsc-c-cuda-k-caching-sp-serial
        COMMAND bin/dwarf-cloudsc-c-cuda-k-caching-sp
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ###

    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )
endif()
else()
    ecbuild_info( "Serialbox and/or CUDA not found, disabling CUDA prototype(s)" )
endif()
