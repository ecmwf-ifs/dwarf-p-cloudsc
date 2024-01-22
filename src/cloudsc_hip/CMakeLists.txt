# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_HIP
	DESCRIPTION "Build the HIP version CLOUDSC using Serialbox" DEFAULT ON
	CONDITION Serialbox_FOUND AND HAVE_HIP
)

if( HAVE_CLOUDSC_HIP )

    set(CMAKE_C_COMPILER "${ROCM_PATH}/bin/hipcc")
    set(CMAKE_CXX_COMPILER "${ROCM_PATH}/bin/hipcc")

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

    if (NOT DEFINED CMAKE_HIP_ARCHITECTURES)
      message(WARNING "No HIP architecture is set! ('CMAKE_HIP_ARCHITECTURES' is not defined)")
    else()
      target_compile_options(dwarf-cloudsc-hip-lib PRIVATE --offload-arch=${CMAKE_HIP_ARCHITECTURES})
    endif()
    
    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-hip
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-hip-lib
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-hip-serial
        COMMAND bin/dwarf-cloudsc-hip
        ARGS 1 1000 64
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
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

    if (NOT DEFINED CMAKE_HIP_ARCHITECTURES)
      message(WARNING "No HIP architecture is set! ('CMAKE_HIP_ARCHITECTURES' is not defined)")
    else()
      target_compile_options(dwarf-cloudsc-hip-hoist-lib PRIVATE --offload-arch=${CMAKE_HIP_ARCHITECTURES})
    endif()

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-hip-hoist
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-hip-hoist-lib
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-hip-hoist-serial
        COMMAND bin/dwarf-cloudsc-hip-hoist
        ARGS 1 1000 64
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
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

    if (NOT DEFINED CMAKE_HIP_ARCHITECTURES)
      message(WARNING "No HIP architecture is set! ('CMAKE_HIP_ARCHITECTURES' is not defined)")
    else()
      target_compile_options(dwarf-cloudsc-hip-k-caching-lib PRIVATE --offload-arch=${CMAKE_HIP_ARCHITECTURES})
    endif()

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-hip-k-caching
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-hip-k-caching-lib
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-hip-k-caching-serial
        COMMAND bin/dwarf-cloudsc-hip-k-caching
        ARGS 1 1000 64
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ## 
    
    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )

endif()
