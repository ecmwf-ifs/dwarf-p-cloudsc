# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define integer IDs corresponding to every language choice
## kernel language
set(CUDA_LANG "1")
set(HIP_LANG "2")
set(SYCL_LANG "3")
## offload model
set(ACC_OFFLOAD "1")
set(OMP_OFFLOAD "2")

# Select offload model
if( HAVE_ACC )
    set( INIT_GPU_OFFLOAD "ACC" )
elseif( HAVE_OMP_TARGET )
    set( INIT_GPU_OFFLOAD "OMP" )
else()
    set( INIT_GPU_OFFLOAD "UNDEFINED" )
endif()

set( CLOUDSC_MIX_GPU_OFFLOAD ${INIT_GPU_OFFLOAD} CACHE STRING "The programming model to use for data offload in mix variants (ACC or OMP)" )

if ( NOT ${CLOUDSC_MIX_GPU_OFFLOAD} MATCHES "ACC|OMP" )
  message(FATAL_ERROR "CLOUDSC_MIX_GPU_OFFLOAD: '${CLOUDSC_MIX_GPU_OFFLOAD}' is not a valid option! (Allowed: 'ACC' and 'OMP')")
endif()
set(GPU_OFFLOAD ${${CLOUDSC_MIX_GPU_OFFLOAD}_OFFLOAD})

# Select kernel language
if( HAVE_CUDA )
    set( INIT_GPU_LANG "CUDA" )
elseif( HAVE_HIP )
    set( INIT_GPU_LANG "HIP" )
elseif( HAVE_SYCL )
    set( INIT_GPU_LANG "SYCL" )
else()
    set( INIT_GPU_LANG "UNDEFINED" )
endif()

set( CLOUDSC_MIX_GPU_LANG ${INIT_GPU_LANG} CACHE STRING "The kernel language to use in mix variants (CUDA, HIP or SYCL)" )

if ( NOT ${CLOUDSC_MIX_GPU_LANG} MATCHES "CUDA|HIP|SYCL" )
    message(FATAL_ERROR "CLOUDSC_MIX_GPU_LANG: '${CLOUDSC_MIX_GPU_LANG}' is not a valid option! (Allowed: 'CUDA', 'HIP' and 'SYCL')")
endif()
set(GPU_LANG ${${CLOUDSC_MIX_GPU_LANG}_LANG})

if (CLOUDSC_MIX_GPU_LANG STREQUAL "CUDA")
    enable_language(CUDA)
    set(CMAKE_CUDA_FLAGS "-O3 -ffast-math")
elseif (CLOUDSC_MIX_GPU_LANG STREQUAL "HIP")
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
elseif (CLOUDSC_MIX_GPU_LANG STREQUAL "SYCL")
    enable_language(CXX)
    set(CMAKE_CXX_STANDARD 17)
endif()

ecbuild_add_library(
    TARGET dwarf-cloudsc-gpu-lib-${prec}
    INSTALL_HEADERS LISTED
    TYPE STATIC
    SOURCES
        yoecldp_c.h
        cloudsc_c_k_caching.cpp
        cloudsc_c_k_caching_mod.F90
    PUBLIC_INCLUDES
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>
    PUBLIC_LIBS
        $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
        $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
        $<${HAVE_CUDA}:CUDA::cudart>
        cloudsc-common-${prec}-lib
)

    set_target_properties( dwarf-cloudsc-gpu-lib-${prec} 
                          PROPERTIES
                          Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module-${CLOUDSC_MIX_GPU_LANG}"
    )


# Provide definitions to target
foreach(_def CUDA_LANG HIP_LANG SYCL_LANG ACC_OFFLOAD OMP_OFFLOAD GPU_LANG GPU_OFFLOAD)
  target_compile_definitions(dwarf-cloudsc-gpu-lib-${prec} PUBLIC ${_def}=${${_def}})
endforeach()

# CUDA specific
if (CLOUDSC_MIX_GPU_LANG STREQUAL "CUDA")
  set_source_files_properties(cloudsc_c_k_caching.cpp PROPERTIES LANGUAGE CUDA)
  target_compile_options( dwarf-cloudsc-gpu-lib-${prec}
      PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:SHELL:${CLOUDSC_CUDA_OPT_FLAGS} ${CLOUDSC_CUDA_FLAGS}>
  )
  set_target_properties( dwarf-cloudsc-gpu-lib-${prec} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
# HIP specific
elseif (CLOUDSC_MIX_GPU_LANG STREQUAL "HIP")
  set_source_files_properties(cloudsc_c_k_caching.cpp PROPERTIES LANGUAGE HIP)
  target_compile_options(dwarf-cloudsc-gpu-lib-${prec} PRIVATE --offload-arch=gfx90a)
# SYCL specific
elseif (CLOUDSC_MIX_GPU_LANG STREQUAL "SYCL")

else ()

endif()

ecbuild_add_executable(
    TARGET  dwarf-cloudsc-mix-${prec}
    SOURCES dwarf_cloudsc_gpu.F90 cloudsc_driver_mod.F90
    LIBS dwarf-cloudsc-gpu-lib-${prec}
)

set_target_properties(dwarf-cloudsc-mix-${prec} PROPERTIES LINKER_LANGUAGE Fortran)
set_target_properties(dwarf-cloudsc-mix-${prec} 
                      PROPERTIES
                      Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module-${CLOUDSC_MIX_GPU_LANG}"
                      )

# CUDA specific
if (CLOUDSC_MIX_GPU_LANG STREQUAL "CUDA")
  set_target_properties( dwarf-cloudsc-mix-${prec} PROPERTIES CUDA_SEPARABLE_COMPILATION ON)
# HIP specific 
elseif (CLOUDSC_MIX_GPU_LANG STREQUAL "HIP")
# SYCL specific
elseif (CLOUDSC_MIX_GPU_LANG STREQUAL "SYCL")
else ()
endif()

ecbuild_add_test(
    TARGET dwarf-cloudsc-mix-${prec}-serial
    COMMAND dwarf-cloudsc-mix-${prec}
    ARGS 1 1000 128
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OMP 1
)

ecbuild_add_test(
    TARGET dwarf-cloudsc-mix-${prec}-mpi
    COMMAND dwarf-cloudsc-mix-${prec}
    ARGS 1 1000 128
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    MPI 2
    OMP 1
    CONDITION HAVE_MPI
)
