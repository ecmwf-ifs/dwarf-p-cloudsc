# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_CROSS
	DESCRIPTION "Build the HIP version CLOUDSC using Serialbox" DEFAULT OFF
	CONDITION Serialbox_FOUND
)

if( HAVE_CLOUDSC_CROSS )
   
    # GPU data offload (default: ACC = 1)
    ## ACC: 1
    ## OMP: 2    
    if (CLOUDSC_GPU_OFFLOAD STREQUAL "ACC")
	set(GPU_OFFLOAD "1")
	# set(OpenACC_Fortran_FLAGS "${OpenACC_Fortran_FLAGS} -gpu=nordc" CACHE STRING "" )
    elseif (CLOUDSC_GPU_OFFLOAD STREQUAL "OMP")
	set(GPU_OFFLOAD "2")
    else ()
	set(GPU_OFFLOAD "1")
    endif()

    # GPU kernel language (default: CUDA = 1)
    ## CUDA: 1
    ## HIP : 2
    ## SYCL: 3
    ##########################################################
    if (CLOUDSC_GPU_LANG STREQUAL "CUDA")
        set(GPU_LANG "1")
	enable_language(CUDA)
	set(CMAKE_CUDA_FLAGS "-O3 -ffast-math")
    ##########################################################
    elseif (CLOUDSC_GPU_LANG STREQUAL "HIP")
	set(GPU_LANG "2")
	if(NOT DEFINED ROCM_PATH)
          if(DEFINED ENV{ROCM_PATH})
            set(ROCM_PATH $ENV{ROCM_PATH} CACHE PATH "Path to which ROCM has been installed")
          else()
            set(ROCM_PATH "/opt/rocm" CACHE PATH "Path to which ROCM has been installed")
          endif()
        endif()
	enable_language(HIP)
        find_package(hip REQUIRED)
	set(CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} -O3 -ffast-math")
    ##########################################################
    elseif (CLOUDSC_GPU_LANG STREQUAL "SYCL")
        set(GPU_LANG "3")
	enable_language(CXX)
        set (CMAKE_CXX_STANDARD 17)
	set(CMAKE_CXX_COMPILER "/home/nams/opt/OpenSYCL/bin/syclcc")
	set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -I/home/nams/opt/OpenSYCL/include --hipsycl-targets=cuda-nvcxx:cc80 -L/home/nams/opt/OpenSYCL/lib")
    ##########################################################
    else ()
        set(GPU_LANG "1")
    endif()
    
    ecbuild_add_library(
        TARGET dwarf-cloudsc-gpu-lib
        INSTALL_HEADERS LISTED
        SOURCES
            load_state.h
    	    load_state.cpp
    	    yoecldp_c.h
    	    cloudsc_c_k_caching.h
      	    cloudsc_c_k_caching.cpp
            cloudsc_c_k_caching_mod.F90
     	    # cloudsc_driver_mod.F90
    	    PUBLIC_INCLUDES
    	    	$<INSTALL_INTERFACE:include>
            	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
    	    PUBLIC_LIBS
    	        cloudsc-common-lib
		# hip::device
    		Serialbox::Serialbox_C
    		# $<${HAVE_OMP}:OpenMP::OpenMP_C>
    		# OpenMP::OpenMP_Fortran
    )

    target_compile_definitions(dwarf-cloudsc-gpu-lib PUBLIC CUDA_LANG=1)
    target_compile_definitions(dwarf-cloudsc-gpu-lib PUBLIC HIP_LANG=2)
    target_compile_definitions(dwarf-cloudsc-gpu-lib PUBLIC SYCL_LANG=3)
    target_compile_definitions(dwarf-cloudsc-gpu-lib PUBLIC ACC_OFFLOAD=1)
    target_compile_definitions(dwarf-cloudsc-gpu-lib PUBLIC OMP_OFFLOAD=2)
    target_compile_definitions(dwarf-cloudsc-gpu-lib PUBLIC GPU_OFFLOAD=${GPU_OFFLOAD})
    target_compile_definitions(dwarf-cloudsc-gpu-lib PUBLIC GPU_LANG=${GPU_LANG})

    target_include_directories(dwarf-cloudsc-gpu-lib PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)

    # CUDA specific
    ##########################################################
    if (CLOUDSC_GPU_LANG STREQUAL "CUDA")
      set_source_files_properties(cloudsc_c_k_caching.cpp load_state.cpp PROPERTIES LANGUAGE CUDA)
      cloudsc_add_compile_options(
         SOURCES cloudsc_c_k_caching_mod.F90 dwarf_cloudsc_gpu.F90 cloudsc_driver_mod.F90
         FLAGS   "-cuda") # -gpu=maxrregcount ...
      target_include_directories(
        dwarf-cloudsc-gpu-lib
        PUBLIC
        ${CMAKE_CUDA_TOOLKIT_INCLUDE_DIRECTORIES}
      )
      if (NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
        target_compile_options(dwarf-cloudsc-gpu-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>>)
      else()
        target_compile_options(dwarf-cloudsc-gpu-lib PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:
              -gencode arch=compute_${CMAKE_CUDA_ARCHITECTURES},code=sm_${CMAKE_CUDA_ARCHITECTURES}>)
      endif()
      set_target_properties( dwarf-cloudsc-gpu-lib PROPERTIES CUDA_SEPARABLE_COMPILATION ON) 
    # HIP specific
    ##########################################################
    elseif (CLOUDSC_GPU_LANG STREQUAL "HIP")
      set_source_files_properties(cloudsc_c_k_caching.cpp load_state.cpp PROPERTIES LANGUAGE HIP) 
    # SYCL specific
    ##########################################################
    elseif (CLOUDSC_GPU_LANG STREQUAL "SYCL")
    
    else ()
    
    endif()
    
    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-cross
	SOURCES dwarf_cloudsc_gpu.F90 cloudsc_driver_mod.F90
	LIBS dwarf-cloudsc-gpu-lib
    )

    # target_compile_definitions(dwarf-cloudsc-cross PUBLIC CUDA_LANG=1)
    # target_compile_definitions(dwarf-cloudsc-cross PUBLIC HIP_LANG=2)
    # target_compile_definitions(dwarf-cloudsc-cross PUBLIC SYCL_LANG=3)
    # target_compile_definitions(dwarf-cloudsc-cross PUBLIC ACC_OFFLOAD=1)
    # target_compile_definitions(dwarf-cloudsc-cross PUBLIC OMP_OFFLOAD=2)
    # target_compile_definitions(dwarf-cloudsc-cross PUBLIC GPU_OFFLOAD=${GPU_OFFLOAD})
    # target_compile_definitions(dwarf-cloudsc-cross PUBLIC GPU_LANG=${GPU_LANG})

    set_target_properties(dwarf-cloudsc-cross PROPERTIES LINKER_LANGUAGE Fortran)

    # CUDA specific
    ##########################################################
    if (CLOUDSC_GPU_LANG STREQUAL "CUDA")
      set_target_properties( dwarf-cloudsc-cross PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
      target_link_libraries(dwarf-cloudsc-cross dwarf-cloudsc-gpu-lib cudadevrt cudart)
      target_link_options(dwarf-cloudsc-cross PRIVATE "-cuda")
    # HIP specific 
    ##########################################################
    elseif (CLOUDSC_GPU_LANG STREQUAL "HIP")
    # SYCL specific
    ##########################################################
    elseif (CLOUDSC_GPU_LANG STREQUAL "SYCL")
    ##########################################################
    else ()
    ##########################################################
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )

endif()
