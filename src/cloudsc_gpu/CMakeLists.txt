# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

### Macro to create executable and accompanying tests

function( cloudsc_add_gpu_exe )

    set( options "" )
    set( single_value_args TARGET COMMON_LIB DRIVER )
    set( multi_value_args DEFINITIONS TEST_ENV KERNEL LIBS )

    cmake_parse_arguments( _CL "${options}" "${single_value_args}" "${multi_value_args}" ${ARGN} )

    if( _CL_UNPARSED_ARGUMENTS )
        ecbuild_critical( "Unknown keywords given to cloudsc_add_gpu_exe(): \"${_CL_UNPARSED_ARGUMENTS}\"")
    endif()

    ecbuild_add_executable(
        TARGET ${_CL_TARGET}-${prec}
        DEFINITIONS ${_CL_DEFINITIONS}
        SOURCES
            dwarf_cloudsc_gpu.F90
            ${_CL_DRIVER}
            ${_CL_KERNEL}
        LIBS
            cloudsc-common-${prec}-${_CL_COMMON_LIB}
            ${_CL_LIBS}
    )

    set_target_properties( ${_CL_TARGET}-${prec}
        PROPERTIES
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module-${_CL_TARGET}"
    )

    ecbuild_add_test(
        TARGET ${_CL_TARGET}-${prec}-serial
        COMMAND ${_CL_TARGET}-${prec}
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        OMP 1
        ENVIRONMENT ${_CL_TEST_ENV}
    )
    ecbuild_add_test(
        TARGET ${_CL_TARGET}-${prec}-omp
        COMMAND ${_CL_TARGET}-${prec}
        ARGS 4 1000 128
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        OMP 4
        ENVIRONMENT ${_CL_TEST_ENV}
        CONDITION HAVE_OMP
    )
    ecbuild_add_test(
        TARGET ${_CL_TARGET}-${prec}-mpi
        COMMAND ${_CL_TARGET}-${prec}
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        MPI 2
        OMP 1
        ENVIRONMENT ${_CL_TEST_ENV}
        CONDITION HAVE_MPI
    )
    ecbuild_add_test(
        TARGET ${_CL_TARGET}-${prec}-mpi-omp
        COMMAND ${_CL_TARGET}-${prec}
        ARGS 2 1000 128
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        MPI 2
        OMP 2
        ENVIRONMENT ${_CL_TEST_ENV}
        CONDITION HAVE_OMP AND HAVE_MPI
    )
endfunction()

set( GPU_VARIANTS scc scc-stack scc-hoist scc-k-caching )
set( CUDA_VARIANTS scc-cuf scc-cuf-k-caching )

### OpenACC GPU variants

if( HAVE_ACC )

    foreach( variant IN LISTS GPU_VARIANTS )
        string( REPLACE "-" "_" variant_ ${variant} )
        string( TOUPPER ${variant_} VARIANT_ )

        if( "${variant}" STREQUAL "scc" )
            set( test_env "NVCOMPILER_ACC_CUDA_HEAPSIZE=64M" )
        else()
            set( test_env )
        endif()

        cloudsc_add_gpu_exe(
            TARGET          dwarf-cloudsc-gpu-${variant}
            DEFINITIONS     CLOUDSC_GPU_${VARIANT_}
            DRIVER          cloudsc_driver_gpu_${variant_}_mod.F90
            KERNEL          cloudsc_gpu_${variant_}_mod.F90
            COMMON_LIB      acc
            TEST_ENV        ${test_env}
        )
    endforeach()

    if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
        target_compile_options(
            dwarf-cloudsc-gpu-scc-stack-${prec}
            PRIVATE "-fcray-pointer"
        )
    elseif( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
        target_compile_options(
            dwarf-cloudsc-gpu-scc-stack-${prec}
            PRIVATE "-Mcray=pointer"
        )
    endif()

endif()

### OpenMP target GPU variants

if( HAVE_OMP AND HAVE_OMP_TARGET )
    set( OMP_DEFINITIONS "CLOUDSC_GPU_OMP_SCC" )
    if( HAVE_OMP_TARGET_LOOP_CONSTRUCT )
        list( PREPEND OMP_DEFINITIONS "HAVE_OMP_TARGET_LOOP_CONSTRUCT" )
    endif()
    if( HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL )
        list( PREPEND OMP_DEFINITIONS "HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL" )
    endif()
    if( HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD )
        list( PREPEND OMP_DEFINITIONS "HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD" )
    endif()

    foreach( variant IN LISTS GPU_VARIANTS )
        string( REPLACE "-" "_" variant_ ${variant} )
        string( TOUPPER ${variant_} VARIANT_ )

        if( "${variant}" STREQUAL "scc" )
            set( test_env "NVCOMPILER_ACC_CUDA_HEAPSIZE=64M" )
        else()
            set( test_env )
        endif()

        cloudsc_add_gpu_exe(
            TARGET      dwarf-cloudsc-gpu-omp-${variant}
            DEFINITIONS "${OMP_DEFINITIONS}_${VARIANT_}"
            DRIVER      cloudsc_driver_gpu_omp_${variant_}_mod.F90
            KERNEL      cloudsc_gpu_omp_${variant_}_mod.F90
            COMMON_LIB  omp
            TEST_ENV    ${test_env}
        )
    endforeach()

    if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU" )
        target_compile_options(
            dwarf-cloudsc-gpu-omp-scc-stack-${prec}
            PRIVATE "-fcray-pointer"
        )
    elseif( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
        target_compile_options(
            dwarf-cloudsc-gpu-omp-scc-stack-${prec}
            PRIVATE "-Mcray=pointer"
        )
    endif()

endif()

### CUDA Fortran variants

if( HAVE_CUDA )

    foreach( variant IN LISTS CUDA_VARIANTS )
        string( REPLACE "-" "_" variant_ ${variant} )
        string( TOUPPER ${variant_} VARIANT_ )

        if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
            cloudsc_add_compile_options(
                SOURCES
                    cloudsc_gpu_${variant_}_mod.F90
                    cloudsc_driver_gpu_${variant_}_mod.F90
                FLAGS "-cuda -gpu=maxregcount:128"
            )
        endif()

        cloudsc_add_gpu_exe(
            TARGET      dwarf-cloudsc-gpu-${variant}
            DEFINITIONS CLOUDSC_GPU_${VARIANT_}
            DRIVER      cloudsc_driver_gpu_${variant_}_mod.F90
            KERNEL
                cloudsc_gpu_${variant_}_mod.F90
                nlev_mod.F90
            COMMON_LIB  cuda
            LIBS        CUDA::cudart
        )
    endforeach()

endif()

### FIELD_API variants

if( HAVE_FIELD_API )

    cloudsc_add_gpu_exe(
        TARGET      dwarf-cloudsc-gpu-scc-field
        DEFINITIONS CLOUDSC_GPU_SCC_FIELD
        DRIVER      cloudsc_driver_gpu_scc_field_mod.F90
        KERNEL      cloudsc_gpu_scc_mod.F90
        COMMON_LIB  field_api
        TEST_ENV
            NVCOMPILER_ACC_CUDA_HEAPSIZE=64M
            DEV_ALLOC_SIZE=1073741824
            $<${HAVE_CUDA}:CLOUDSC_FIELD_API_PINNED=ON>
    )

    # Importantly, we add the "pinned" flag to the linker command to ensure pinning!
    if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" AND field_api_HAVE_CUDA )
        target_link_options( dwarf-cloudsc-gpu-scc-field-${prec} PUBLIC "-gpu=pinned" )
    endif()
endif()
