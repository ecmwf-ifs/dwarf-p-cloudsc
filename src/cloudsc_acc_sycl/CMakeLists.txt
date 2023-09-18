# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_ACC_SYCL
	DESCRIPTION "Build the CUDA version CLOUDSC using Serialbox" DEFAULT OFF
        CONDITION Serialbox_FOUND
)

if( HAVE_CLOUDSC_ACC_SYCL )

   enable_language(CXX)
   set (CMAKE_CXX_STANDARD 17)

   # find_package(hipSYCL REQUIRED)

   # enable_language(CUDA)
   # set(CMAKE_CUDA_FLAGS "-O3 -ffast-math")

   set(CMAKE_CXX_COMPILER "/home/nams/opt/OpenSYCL/bin/syclcc")
   # --hipsycl-clang=/usr/bin/clang
   # --hipsycl-targets=cuda:sm_80
   #set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} --hipsycl-platform=cuda --hipsycl-targets=cuda-nvcxx --hipsycl-nvcxx=/home/nams/opt/OpenSYCL/bin/syclcc -fcuda-approx-transcendentals -fcuda-flush-denormals-to-zero --no-cuda-version-check --cuda-path=/usr/local/apps/nvidia/22.1/Linux_x86_64/22.1/cuda")

   set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -I/home/nams/opt/OpenSYCL/include --hipsycl-targets=cuda-nvcxx:cc80 -L/home/nams/opt/OpenSYCL/lib")

    ###### SCC-HIP ####

    ecbuild_add_library(
	    TARGET dwarf-cloudsc-sycl-lib
	    INSTALL_HEADERS LISTED
	    SOURCES
	        load_state.h
		load_state.cpp
		yoecldp_c.h
		# cloudsc_c_k_caching.h
		cloudsc_c_k_caching.cpp
		cloudsc_c_k_caching.kernel
		cloudsc_c_k_caching_mod.F90
		# cloudsc_driver_mod.F90
	    PUBLIC_INCLUDES
	    	$<INSTALL_INTERFACE:include>
            	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
	    PUBLIC_LIBS
	    	cloudsc-common-lib
		Serialbox::Serialbox_C
		# $<${HAVE_OMP}:OpenMP::OpenMP_C>
		# OpenMP::OpenMP_Fortran
		# TYPE STATIC
    )

    # cloudsc_add_compile_options(
    #      SOURCES cloudsc_c_k_caching_mod.F90 dwarf_cloudsc_gpu.F90 cloudsc_driver_mod.F90
    #      FLAGS   "-cuda") # -gpu=maxrregcount ... 
    # target_link_options(dwarf-cloudsc-cuda-lib PRIVATE "-cuda --device-link")

    target_include_directories(dwarf-cloudsc-sycl-lib PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
    
    # link_directories(/usr/local/apps/nvidia/22.11/Linux_x86_64/22.11/cuda/lib64)
    # target_link_libraries(dwarf-cloudsc-cuda-lib PUBLIC cudadevrt cudart)
    
    target_include_directories(
      dwarf-cloudsc-sycl-lib
      PUBLIC
      ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
    )
    # if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
    #   target_compile_options(dwarf-cloudsc-sycl-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>>)
    # else()
    #   target_compile_options(dwarf-cloudsc-sycl-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
    #         -gencode arch=compute_${CMAKE_CUDA_ARCHITECTURES},code=sm_${CMAKE_CUDA_ARCHITECTURES}>)
    # endif()

    # set_target_properties( dwarf-cloudsc-sycl-lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-acc-sycl
	SOURCES dwarf_cloudsc_gpu.F90 cloudsc_driver_mod.F90
        LIBS dwarf-cloudsc-sycl-lib
    )

    # set_target_properties( dwarf-cloudsc-acc-sycl PROPERTIES CUDA_SEPARABLE_COMPILATION ON)

    set_target_properties(dwarf-cloudsc-acc-sycl PROPERTIES LINKER_LANGUAGE Fortran)
    # target_link_libraries(dwarf-cloudsc-acc-sycl /home/nams/opt/OpenSYCL/lib/dwarf-cloudsc-sycl-lib.so)# libboost_filesystem libboost_system)
    # target_link_options(dwarf-cloudsc-acc-sycl PRIVATE "-cuda")

    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )

endif()
