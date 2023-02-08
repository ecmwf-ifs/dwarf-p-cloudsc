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


    enable_language(CUDA)
    enable_language(CXX)

    ecbuild_add_library(
        TARGET dwarf-cloudsc-cuda-lib
	INSTALL_HEADERS LISTED
        SOURCES
	    cloudsc/yoecldp_c.h
	    cloudsc/yoethf_c.h
	    cloudsc/yomcst_c.h
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


    target_compile_options(dwarf-cloudsc-cuda-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
    #                   --generate-line-info
    #                   --relocatable-device-code=true
    #		       -g -G
                       -gencode arch=compute_80,code=sm_80
                       >)

    target_include_directories(dwarf-cloudsc-cuda-lib PUBLIC "/usr/local/apps/nvidia/22.1/Linux_x86_64/22.1/cuda/include")
   
    #target_link_options(dwarf-cloudsc-cuda-lib PUBLIC "-Mcuda") 
    set_target_properties( dwarf-cloudsc-cuda-lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
    
    #target_include_directories( dwarf-cloudsc-cuda-lib PUBLIC
    #    ${CMAKE_CURRENT_SOURCE_DIR}
    #    ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    #    ${CUDA_INCLUDE_DIRS})

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-cuda
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-cuda-lib
    )

    target_link_libraries(dwarf-cloudsc-cuda dwarf-cloudsc-cuda-lib)

    ###### SCC-CUDA-HOIST ####
    ecbuild_add_library(
        TARGET dwarf-cloudsc-cuda-hoist-lib
        INSTALL_HEADERS LISTED
        SOURCES
            cloudsc/yoecldp_c.h
            cloudsc/yoethf_c.h
            cloudsc/yomcst_c.h
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


    target_compile_options(dwarf-cloudsc-cuda-hoist-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
    #                   --generate-line-info
    #                   --relocatable-device-code=true
    #                  -g -G
                       -gencode arch=compute_80,code=sm_80
                       >)

    target_include_directories(dwarf-cloudsc-cuda-hoist-lib PUBLIC "/usr/local/apps/nvidia/22.1/Linux_x86_64/22.1/cuda/include")

    #target_link_options(dwarf-cloudsc-cuda-lib PUBLIC "-Mcuda") 
    set_target_properties( dwarf-cloudsc-cuda-hoist-lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

    #target_include_directories( dwarf-cloudsc-cuda-lib PUBLIC
    #    ${CMAKE_CURRENT_SOURCE_DIR}
    #    ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    #    ${CUDA_INCLUDE_DIRS})

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-cuda-hoist
        SOURCES dwarf_cloudsc_hoist.cpp
        LIBS dwarf-cloudsc-cuda-hoist-lib
    )

    target_link_libraries(dwarf-cloudsc-cuda-hoist dwarf-cloudsc-cuda-hoist-lib)
    ###


    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-cuda-serial
        COMMAND bin/dwarf-cloudsc-cuda
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-cuda-omp
        COMMAND bin/dwarf-cloudsc-cuda
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 4
        CONDITION HAVE_OMP
    )

else()
    ecbuild_info( "Serialbox not found, disabling C prototype" )
endif()