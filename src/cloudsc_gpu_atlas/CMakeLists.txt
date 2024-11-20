# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


ecbuild_add_option( FEATURE CLOUDSC_GPU_SCC
    DESCRIPTION "Build optimized GPU version of CLOUDSC based on Atlas + OpenACC" DEFAULT ON
    CONDITION atlas_FOUND AND (Serialbox_FOUND OR HDF5_FOUND)
)


ecbuild_add_option( FEATURE CLOUDSC_GPU_SCC_HOIST
    DESCRIPTION "Build optimized GPU version of CLOUDSC using SCC with hoisted temporary arrays (based on Atlas + OpenACC)" DEFAULT ON
    CONDITION atlas_FOUND AND (Serialbox_FOUND OR HDF5_FOUND)
)

ecbuild_add_option( FEATURE CLOUDSC_GPU_SCC_K_CACHING
    DESCRIPTION "Build optimized GPU version of CLOUDSC using SCC with k-caching of local array in the kernel (based on Atlas + OpenACC)" DEFAULT ON
    CONDITION atlas_FOUND AND (Serialbox_FOUND OR HDF5_FOUND)
)


if( HAVE_CLOUDSC_GPU_SCC )
    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-scc-atlas
        SOURCES
            cloudsc_global_atlas_state_mod.F90
            expand_atlas_mod.F90
            validate_atlas_mod.F90
            ../cloudsc_gpu/cloudsc_gpu_scc_mod.F90
            cloudsc_driver_scc_mod.F90
            dwarf_cloudsc_atlas.F90
        LIBS
            atlas_f
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-gpu-scc-atlas
        COMMAND bin/dwarf-cloudsc-gpu-scc-atlas
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
        ENVIRONMENT "NVCOMPILER_ACC_CUDA_HEAPSIZE=64M"
    )
endif()


if( HAVE_CLOUDSC_GPU_SCC_HOIST )
    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-scc-hoist-atlas
        SOURCES
            cloudsc_global_atlas_state_mod.F90
            expand_atlas_mod.F90
            validate_atlas_mod.F90
            ../cloudsc_gpu/cloudsc_gpu_scc_hoist_mod.F90
            cloudsc_driver_scc_mod.F90
            dwarf_cloudsc_atlas.F90
        LIBS
            atlas_f
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_HOIST
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-gpu-scc-hoist-atlas
        COMMAND bin/dwarf-cloudsc-gpu-scc-hoist-atlas
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
endif()


if( HAVE_CLOUDSC_GPU_SCC_K_CACHING )
    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-scc-k-caching-atlas
        SOURCES
            cloudsc_global_atlas_state_mod.F90
            expand_atlas_mod.F90
            validate_atlas_mod.F90
            ../cloudsc_gpu/cloudsc_gpu_scc_k_caching_mod.F90
            cloudsc_driver_scc_mod.F90
            dwarf_cloudsc_atlas.F90
        LIBS
            atlas_f
            cloudsc-common-lib
            DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_K_CACHING
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-gpu-scc-k-caching-atlas
        COMMAND bin/dwarf-cloudsc-gpu-scc-k-caching-atlas
        ARGS 1 1000 128
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
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


# Importantly, we add the "pinned" flag to the linker command to ensure pinning!
#if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
#    target_link_options( dwarf-cloudsc-gpu-scc-atlas PUBLIC "-gpu=pinned" )
#endif()
