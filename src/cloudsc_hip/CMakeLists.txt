# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_CUDA
    DESCRIPTION "Build the CUDA version CLOUDSC using Serialbox" DEFAULT ON
    CONDITION Serialbox_FOUND
)

if( HAVE_CLOUDSC_CUDA )


    if(NOT DEFINED HIP_PATH)
      if(NOT DEFINED ENV{HIP_PATH})
        set(HIP_PATH "/opt/rocm-5.0.2/hip" CACHE PATH "Path to which HIP has been installed")
      else()
        set(HIP_PATH $ENV{HIP_PATH} CACHE PATH "Path to which HIP has been installed")
      endif()
    endif()

    if(NOT DEFINED ROCM_PATH)
      if(DEFINED ENV{ROCM_PATH})
        set(ROCM_PATH $ENV{ROCM_PATH} CACHE PATH "Path to which ROCM has been installed")
      else()
        set(ROCM_PATH "/opt/rocm-5.0.2" CACHE PATH "Path to which ROCM has been installed")
      endif()
    endif()

    if ( HAVE_OMP )
      message("found OpenMP for HIP!!")
    else()
      # message( FATAL_ERROR "have NOT found OpenMP for HIP!!")
      message("have not found OMP for HIP!")
    endif()

    set(CMAKE_MODULE_PATH "${HIP_PATH}/cmake" ${CMAKE_MODULE_PATH})
    list(APPEND CMAKE_PREFIX_PATH
      "${HIP_PATH}/lib/cmake"
      "${HIP_PATH}/../lib/cmake" # hopefully catches all extra HIP dependencies
    )
    list(APPEND CMAKE_PREFIX_PATH ${HIP_PATH} ${ROCM_PATH})
    
    find_package(HIP QUIET)
    if(HIP_FOUND)
      message(STATUS "Found HIP: " ${HIP_VERSION})
    else()
      message(FATAL_ERROR "Could not find HIP. Ensure that HIP is either installed in /opt/rocm/hip or the variable HIP_PATH is set to point to the right location.")
    endif()
    find_package(hip REQUIRED)

    set(CMAKE_C_COMPILER "/cfs/klemming/root/rocm/opt/rocm-5.3.3/bin/hipcc")
    set(CMAKE_CXX_COMPILER "/cfs/klemming/root/rocm/opt/rocm-5.3.3/bin/hipcc")
    set(CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -O3")
    set(CMAKE_CXX_LINKER "/cfs/klemming/root/rocm/opt/rocm-5.3.3/bin/hipcc")
    set(CMAKE_HIPCXX_FLAGS ${CMAKE_HIPCXX_FLAGS} -DHIP_PLATFORM=amd -DHAVE_HIP -DHIP_FAST_MATH -O3 --offload-arch=gfx90a )

    ###### SCC-HIP ####
    ecbuild_add_library(
        TARGET dwarf-cloudsc-hip-lib
    	INSTALL_HEADERS LISTED
        SOURCES
    	    cloudsc/yoecldp_c.h
    	    cloudsc/load_state.h
            cloudsc/load_state.cpp
            cloudsc/cloudsc_c.h
            cloudsc/cloudsc_c.cpp
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
    	    hip::device
            Serialbox::Serialbox_C
            $<${HAVE_OMP}:OpenMP::OpenMP_C>
    )

    target_include_directories(dwarf-cloudsc-hip-lib PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>)
    target_link_libraries(dwarf-cloudsc-hip-lib PUBLIC hip::device Serialbox::Serialbox_C $<${HAVE_OMP}:OpenMP::OpenMP_C>)

    target_compile_options(dwarf-cloudsc-hip-lib PRIVATE --offload-arch=gfx90a)
    
    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-hip
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-hip-lib
    )

    ## 

    ###### SCC-HOIST-HIP ####
    ecbuild_add_library(
        TARGET dwarf-cloudsc-hip-hoist-lib
    	INSTALL_HEADERS LISTED
        SOURCES
    	    cloudsc/yoecldp_c.h
    	    cloudsc/load_state.h
            cloudsc/load_state.cpp
            cloudsc/cloudsc_c_hoist.h
            cloudsc/cloudsc_c_hoist.cpp
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
    	    hip::device
            Serialbox::Serialbox_C
            $<${HAVE_OMP}:OpenMP::OpenMP_C>
    )

    target_include_directories(dwarf-cloudsc-hip-hoist-lib PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>)
    target_link_libraries(dwarf-cloudsc-hip-hoist-lib PUBLIC hip::device Serialbox::Serialbox_C $<${HAVE_OMP}:OpenMP::OpenMP_C>)

    target_compile_options(dwarf-cloudsc-hip-hoist-lib PRIVATE --offload-arch=gfx90a)
    
    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-hip-hoist
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-hip-hoist-lib
    )

    ## 

    ###### SCC-K-CACHING-HIP ####
    ecbuild_add_library(
        TARGET dwarf-cloudsc-hip-k-caching-lib
    	INSTALL_HEADERS LISTED
        SOURCES
    	    cloudsc/yoecldp_c.h
    	    cloudsc/load_state.h
            cloudsc/load_state.cpp
            cloudsc/cloudsc_c_k_caching.h
            cloudsc/cloudsc_c_k_caching.cpp
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
    	    hip::device
            Serialbox::Serialbox_C
            $<${HAVE_OMP}:OpenMP::OpenMP_C>
    )

    target_include_directories(dwarf-cloudsc-hip-k-caching-lib PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>)
    target_link_libraries(dwarf-cloudsc-hip-k-caching-lib PUBLIC hip::device Serialbox::Serialbox_C $<${HAVE_OMP}:OpenMP::OpenMP_C>)

    target_compile_options(dwarf-cloudsc-hip-k-caching-lib PRIVATE --offload-arch=gfx90a)
    
    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-hip-k-caching
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-hip-k-caching-lib
    )

    ## 

    
    
    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-hip-serial
        COMMAND bin/dwarf-cloudsc-hip
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-hip-omp
        COMMAND bin/dwarf-cloudsc-hip
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 4
        CONDITION HAVE_OMP
    )

#else()
#    ecbuild_info( "Serialbox not found, disabling C prototype" )
endif()
