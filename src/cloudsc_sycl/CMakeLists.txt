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
    CONDITION Serialbox_FOUND
)

if( HAVE_CLOUDSC_SYCL )

    # set( CMAKE_CXX_STANDARD 17 )
    # set( CMAKE_CXX_STANDARD_REQUIRED ON )
    #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} /EHsc")

    #find_package(LLVM REQUIRED CONFIG)
    #list(APPEND CMAKE_MODULE_PATH "/usr/local/apps/intel/2021.4.0/compiler/latest/linux/cmake/SYCL")
    #find_package(IntelDPCPP REQUIRED)
    # set(LLVM_TARGETS_TO_BUILD "X86;NVPTX")
    
    #set(LLVM_ALL_TARGETS
    #    NVPTX
    #    X86
    #   )

    #set(LLVM_TARGETS_TO_BUILD "all"
    #CACHE STRING "Semicolon-separated list of targets to build, or \"all\".")

    #if( LLVM_TARGETS_TO_BUILD STREQUAL "all" )
    #  set( LLVM_TARGETS_TO_BUILD ${LLVM_ALL_TARGETS} )
    #endif()
    
    #set(LLVM_TARGETS_TO_BUILD
    #  ${LLVM_TARGETS_TO_BUILD})
    #list(REMOVE_DUPLICATES LLVM_TARGETS_TO_BUILD)

    set(CMAKE_CXX_COMPILER "/home/nams/opt/dpcpp/bin/clang++")
    set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/local/apps/intel/2021.4.0/compiler/2021.4.0/linux/compiler/include -O3 -L/home/nams/opt/dpcpp/lib -fopenmp -fsycl-early-optimizations -fsycl -fsycl-targets=nvptx64-nvidia-cuda -Xcuda-ptxas --maxrregcount=128 -Xsycl-target-backend --cuda-gpu-arch=sm_80")
    
    
    ecbuild_add_library(
        TARGET dwarf-cloudsc-cuda-sycl-lib
        INSTALL_HEADERS LISTED
        SOURCES
            cloudsc/yoecldp_c.h
	    cloudsc/load_state.h
            cloudsc/load_state.cpp
	    #cloudsc/cloudsc_c.h
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
	    sycl
            Serialbox::Serialbox_C
    	    $<${HAVE_OMP}:OpenMP::OpenMP_C>
    )


    # target_include_directories(dwarf-cloudsc-cuda-sycl-lib PRIVATE /usr/local/apps/intel/2021.4.0/compiler/2021.4.0/linux/compiler/include)
    # target_compile_options(dwarf-cloudsc-sycl-lib -fsycl-targets=spir64_gen-unknown-unknown-sycldevice) 
    # target_compile_options(dwarf-cloudsc-sycl-lib PUBLIC -fsycl -fsycl-targets=nvptx64-nvidia-cuda -fno-sycl-libspirv) # -fsycl-targets=nvptx64-*-*-sycldevice,spir64-*-*-sycldevice -fno-sycl-libspirv --cuda-path=/usr/local/apps/nvidia/22.1/Linux_x86_64/22.1/cuda/11.5/nvvm/libdevice -Xcuda-ptxas --verbose -Xsycl-target-backend --cuda-gpu-arch=sm_80)#-fsycl-targets=nvptx64-nvidia-cuda -fno-sycl-libspirv -Xcuda-ptxas --maxrregcount=128 -Xcuda-ptxas --verbose -Xsycl-target-backend --cuda-gpu-arch=sm_80) #nvidia_gpu_sm_80) #spir64_gen-unknown-unknown-sycldevice) # -fsycl-targets=nvptx64-nvidia-cuda-sycldevice) # -lpthreads)# -fsycl-targets=nvptx64-nvidia-cuda-sycldevice)
    #target_link_libraries(dwarf-cloudsc-sycl-lib sycl)
    # target_link_options(dwarf-cloudsc-sycl-lib PUBLIC -fsycl -fsycl-targets=nvptx64-nvidia-cuda -fno-sycl-libspirv) # -fsycl-targets=nvptx64-*-*-sycldevice,spir64-*-*-sycldevice -fno-sycl-libspirv --cuda-path=/usr/local/apps/nvidia/22.1/Linux_x86_64/22.1/cuda/11.5/nvvm/libdevice -Xcuda-ptxas --verbose -Xsycl-target-backend --cuda-gpu-arch=sm_80) #-fsycl-targets=nvptx64-nvidia-cuda -fno-sycl-libspirv -Xcuda-ptxas --maxrregcount=128 -Xcuda-ptxas --verbose -Xsycl-target-backend --cuda-gpu-arch=sm_80) #nvidia_gpu_sm_80) #spir64_gen-unknown-unknown-sycldevice) #-fsycl-targets=nvptx64-nvidia-cuda-sycldevice) # -lpthreads)# -fsycl-targets=nvptx64-nvidia-cuda-sycldevice)

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-cuda-sycl
        SOURCES dwarf_cloudsc.cpp
        LIBS dwarf-cloudsc-cuda-sycl-lib
    )

    #target_compile_options(dwarf-cloudsc-sycl PUBLIC -fsycl)# -fsycl-targets=nvptx64-nvidia-cuda-sycldevice)
    # target_link_libraries(dwarf-cloudsc-sycl sycl)
    #target_link_options(dwarf-cloudsc-sycl PUBLIC -fsycl)# -fsycl-targets=nvptx64-nvidia-cuda-sycldevice)

    # Create symlink for the input data
    #execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
    #    ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )

    #ecbuild_add_test(
    #    TARGET dwarf-cloudsc-sycl-serial
    #    COMMAND bin/dwarf-cloudsc-c
    #    ARGS 1 100 16
    #    WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
    #    OMP 1
    #)

    # ecbuild_add_executable(
    #	    TARGET hello
    #	    SOURCES hello.cpp
    # 	    #INCLUDES 
    #	    #    /usr/local/apps/intel/2021.4.0/compiler/latest/linux/include/sycl
    #	    #		/usr/local/apps/intel/2021.4.0/compiler/latest/linux/include
    #)

    #target_compile_options(hello PUBLIC -fsycl)# -fsycl-targets=nvptx64-nvidia-cuda-sycldevice)
    #target_link_libraries(hello sycl)
    #target_link_options(hello PUBLIC -fsycl)# -fsycl-targets=nvptx64-nvidia-cuda-sycldevice)

else()
	ecbuild_info( "Serialbox not found, disabling SYCL version" )
endif()
