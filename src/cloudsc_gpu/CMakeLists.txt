# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define the naive GPU dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_GPU_KERNELS
    DESCRIPTION "Build naive GPU version of CLOUDSC using `!$acc kernels`" DEFAULT OFF
    CONDITION Serialbox_FOUND OR HDF5_FOUND
)

# Define the CLAW-based GPU dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_GPU_CLAW
    DESCRIPTION "Build optimized GPU version of CLOUDSC derived from CLAW using OpenACC" DEFAULT OFF
    CONDITION Serialbox_FOUND OR HDF5_FOUND
)

ecbuild_add_option( FEATURE CLOUDSC_GPU_SCC
    DESCRIPTION "Build optimized GPU version of CLOUDSC using SCC layout and OpenACC" DEFAULT OFF
    CONDITION Serialbox_FOUND OR HDF5_FOUND
)

ecbuild_add_option( FEATURE CLOUDSC_GPU_SCC_CUF
    DESCRIPTION "Build optimized GPU version of CLOUDSC using SCC layout in CUDA" DEFAULT OFF
    CONDITION HAVE_CUDA AND ( Serialbox_FOUND OR HDF5_FOUND )
)

ecbuild_add_option( FEATURE CLOUDSC_GPU_SCC_CUF_K_CACHING
    DESCRIPTION "Build (further) optimized GPU version of CLOUDSC using SCC layout in CUDA" DEFAULT OFF
    CONDITION HAVE_CUDA AND ( Serialbox_FOUND OR HDF5_FOUND )
)

ecbuild_add_option( FEATURE CLOUDSC_GPU_SCC_HOIST
    DESCRIPTION "Build optimized GPU version of CLOUDSC using SCC with hoisted temporary arrays and OpenACC" DEFAULT OFF
    CONDITION Serialbox_FOUND OR HDF5_FOUND
)

ecbuild_add_option( FEATURE CLOUDSC_GPU_SCC_K_CACHING
    DESCRIPTION "Build (further) optimized GPU version of CLOUDSC using SCC layout with OpenACC" DEFAULT OFF
    CONDITION Serialbox_FOUND OR HDF5_FOUND
)

ecbuild_add_option( FEATURE CLOUDSC_GPU_OMP_SCC_HOIST
    DESCRIPTION "Build optimized GPU version of CLOUDSC using SCC with hoisted temporary arrays and OpenMP offload" DEFAULT OFF
    CONDITION Serialbox_FOUND OR HDF5_FOUND
)

ecbuild_add_option( FEATURE CLOUDSC_GPU_SCC_FIELD
    DESCRIPTION "Build optimized GPU version of CLOUDSC using SCC with FIELD API" DEFAULT ON
    CONDITION HAVE_FIELD_API AND ( Serialbox_FOUND OR HDF5_FOUND )
)


if( HAVE_CLOUDSC_GPU_CLAW )
    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-claw
        SOURCES
            dwarf_cloudsc_gpu.F90
            cloudsc_driver_gpu_claw_mod.F90
            cloudsc.claw.gpu.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_CLAW
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-gpu-claw-serial
        COMMAND bin/dwarf-cloudsc-gpu-claw
        ARGS 1 1280 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
        ENABLED OFF  # CLAW variant is currently broken
    )
endif()


if( HAVE_CLOUDSC_GPU_SCC )
    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-scc
        SOURCES
            dwarf_cloudsc_gpu.F90
            cloudsc_driver_gpu_scc_mod.F90
            cloudsc_gpu_scc_mod.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-gpu-scc-serial
        COMMAND bin/dwarf-cloudsc-gpu-scc
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
        ENVIRONMENT "NVCOMPILER_ACC_CUDA_HEAPSIZE=64M"
    )
endif()


if( HAVE_CLOUDSC_GPU_SCC_HOIST )
    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-scc-hoist
        SOURCES
            dwarf_cloudsc_gpu.F90
            cloudsc_driver_gpu_scc_hoist_mod.F90
            cloudsc_gpu_scc_hoist_mod.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_HOIST
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-gpu-scc-hoist-serial
        COMMAND bin/dwarf-cloudsc-gpu-scc-hoist
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
endif()

if( HAVE_CLOUDSC_GPU_SCC_K_CACHING )
    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-scc-k-caching
        SOURCES
            dwarf_cloudsc_gpu.F90
            cloudsc_driver_gpu_scc_k_caching_mod.F90
            cloudsc_gpu_scc_k_caching_mod.F90
        LIBS
            cloudsc-common-lib
	DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_K_CACHING
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-gpu-scc-k-caching-serial
        COMMAND bin/dwarf-cloudsc-gpu-scc-k-caching
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
endif()

if( HAVE_CLOUDSC_GPU_OMP_SCC_HOIST )
    list( APPEND CLOUDSC_GPU_OMP_SCC_HOIST_DEFINITIONS CLOUDSC_GPU_OMP_SCC_HOIST )
    if( HAVE_OMP_TARGET_LOOP_CONSTRUCT )
        list( APPEND CLOUDSC_GPU_OMP_SCC_HOIST_DEFINITIONS HAVE_OMP_TARGET_LOOP_CONSTRUCT )
        if( HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL )
            list( APPEND CLOUDSC_GPU_OMP_SCC_HOIST_DEFINITIONS HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL )
        endif()
        if( HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD )
            list( APPEND CLOUDSC_GPU_OMP_SCC_HOIST_DEFINITIONS HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD )
        endif()
    endif()

    if( NOT (HAVE_OMP_TARGET_TEAMS_DISTRIBUTE OR HAVE_OMP_TARGET_LOOP_CONSTRUCT) )
        ecbuild_warn( "No support for OpenMP target offload, binary will be CPU-only" )
    endif()

    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-omp-scc-hoist
        SOURCES
            dwarf_cloudsc_gpu.F90
            cloudsc_driver_gpu_omp_scc_hoist_mod.F90
            cloudsc_gpu_omp_scc_hoist_mod.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} ${CLOUDSC_GPU_OMP_SCC_HOIST_DEFINITIONS}
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-gpu-omp-scc-hoist-serial
        COMMAND bin/dwarf-cloudsc-gpu-omp-scc-hoist
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
endif()


if( HAVE_CLOUDSC_GPU_SCC_CUF )
    # Compile CUDA fortran files with -cuda.
    cloudsc_add_compile_options(
        SOURCES cloudsc_gpu_scc_cuf_mod.F90  cloudsc_driver_gpu_scc_cuf_mod.F90
        FLAGS   "-cuda -gpu=maxregcount:128")

    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-scc-cuf
        SOURCES
            dwarf_cloudsc_gpu.F90
	    cloudsc_driver_gpu_scc_cuf_mod.F90
	    cloudsc_gpu_scc_cuf_mod.F90
	    nlev_mod.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_CUF
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-gpu-scc-cuf-serial
        COMMAND bin/dwarf-cloudsc-gpu-scc-cuf
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
endif()

if ( HAVE_CLOUDSC_GPU_SCC_CUF_K_CACHING )
    # NEW CUF with k-caching!!!!
    # Compile CUDA fortran files with -cuda.
    cloudsc_add_compile_options(
        SOURCES cloudsc_gpu_scc_cuf_k_caching_mod.F90  cloudsc_driver_gpu_scc_cuf_k_caching_mod.F90
        FLAGS   "-cuda -gpu=maxregcount:128")

    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-scc-cuf-k-caching
        SOURCES
            dwarf_cloudsc_gpu.F90
            cloudsc_driver_gpu_scc_cuf_k_caching_mod.F90
            cloudsc_gpu_scc_cuf_k_caching_mod.F90
            nlev_mod.F90
        LIBS
            cloudsc-common-lib
	    DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_CUF_K_CACHING
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-gpu-scc-cuf-k-caching-serial
        COMMAND bin/dwarf-cloudsc-gpu-scc-cuf-k-caching
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )

endif()


if( HAVE_CLOUDSC_GPU_SCC_FIELD )
    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-scc-field
        SOURCES
            dwarf_cloudsc_gpu.F90
            cloudsc_driver_gpu_scc_field_mod.F90
            cloudsc_gpu_scc_mod.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_FIELD
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-gpu-scc-field
        COMMAND bin/dwarf-cloudsc-gpu-scc-field
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
        ENVIRONMENT "NVCOMPILER_ACC_CUDA_HEAPSIZE=64M"
    )

    # Importantly, we add the "pinned" flag to the linker command to ensure pinning!
    if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
        target_link_options( dwarf-cloudsc-gpu-scc-field PUBLIC "-gpu=pinned" )
    endif()
endif()


# Create symlink for the input data
if( HAVE_SERIALBOX )
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )
endif()

if( HAVE_HDF5 )
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../input.h5 )
  execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../reference.h5 )
endif()
