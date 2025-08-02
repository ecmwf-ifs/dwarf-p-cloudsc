# (C) Copyright 2025- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


ecbuild_add_executable(
    TARGET dwarf-cloudsc-gpu-scc-atlas-${prec}
    SOURCES
        cloudsc_global_atlas_state_mod.F90
        ../cloudsc_fortran_atlas/expand_atlas_mod.F90
        ../cloudsc_fortran_atlas/validate_atlas_mod.F90
        ../cloudsc_gpu/cloudsc_gpu_scc_mod.F90
        cloudsc_driver_scc_mod.F90
        dwarf_cloudsc_atlas.F90
    LIBS
        atlas_f
        pluto_f
        cloudsc-common-${prec}-atlas
    DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC
)

set_target_properties( dwarf-cloudsc-gpu-scc-atlas-${prec}
        PROPERTIES
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module-dwarf-cloudsc-gpu-scc-atlas"
)

ecbuild_add_test(
    TARGET dwarf-cloudsc-gpu-scc-atlas-${prec}
    COMMAND bin/dwarf-cloudsc-gpu-scc-atlas-${prec}
    ARGS 1 1000 128
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OMP 1
    ENVIRONMENT "NVCOMPILER_ACC_CUDA_HEAPSIZE=64M" CLOUDSC_ATLAS_MULTIFIELD=0
)

ecbuild_add_test(
    TARGET dwarf-cloudsc-gpu-scc-multifield-atlas-${prec}
    COMMAND bin/dwarf-cloudsc-gpu-scc-atlas-${prec}
    ARGS 1 1000 128
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OMP 1
    ENVIRONMENT "NVCOMPILER_ACC_CUDA_HEAPSIZE=64M" CLOUDSC_ATLAS_MULTIFIELD=1
)

ecbuild_add_executable(
    TARGET dwarf-cloudsc-gpu-scc-hoist-atlas-${prec}
    SOURCES
        cloudsc_global_atlas_state_mod.F90
        ../cloudsc_fortran_atlas/expand_atlas_mod.F90
        ../cloudsc_fortran_atlas/validate_atlas_mod.F90
        ../cloudsc_gpu/cloudsc_gpu_scc_hoist_mod.F90
        cloudsc_driver_scc_mod.F90
        dwarf_cloudsc_atlas.F90
    LIBS
        atlas_f
        pluto_f
        cloudsc-common-${prec}-atlas
    DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_HOIST
)

set_target_properties( dwarf-cloudsc-gpu-scc-hoist-atlas-${prec}
        PROPERTIES
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module-dwarf-cloudsc-gpu-scc-hoist-atlas"
)

ecbuild_add_test(
    TARGET dwarf-cloudsc-gpu-scc-hoist-atlas-${prec}
    COMMAND bin/dwarf-cloudsc-gpu-scc-hoist-atlas-${prec}
    ARGS 1 1000 128
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OMP 1
    ENVIRONMENT HOIST_POOL=0 CLOUDSC_ATLAS_MULTIFIELD=0
)

ecbuild_add_test(
    TARGET dwarf-cloudsc-gpu-scc-hoist-multifield-atlas-${prec}
    COMMAND bin/dwarf-cloudsc-gpu-scc-hoist-atlas-${prec}
    ARGS 1 1000 128
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OMP 1
    ENVIRONMENT HOIST_POOL=0 CLOUDSC_ATLAS_MULTIFIELD=1
)

ecbuild_add_test(
    TARGET dwarf-cloudsc-gpu-scc-hoist_pool-atlas-${prec}
    COMMAND bin/dwarf-cloudsc-gpu-scc-hoist-atlas-${prec}
    ARGS 1 1000 128
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OMP 1
    ENVIRONMENT HOIST_POOL=1 CLOUDSC_ATLAS_MULTIFIELD=0
)

ecbuild_add_test(
    TARGET dwarf-cloudsc-gpu-scc-hoist_pool-multifield-atlas-${prec}
    COMMAND bin/dwarf-cloudsc-gpu-scc-hoist-atlas-${prec}
    ARGS 1 1000 128
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OMP 1
    ENVIRONMENT HOIST_POOL=1 CLOUDSC_ATLAS_MULTIFIELD=1
)

ecbuild_add_executable(
    TARGET dwarf-cloudsc-gpu-scc-k-caching-atlas-${prec}
    SOURCES
        cloudsc_global_atlas_state_mod.F90
        ../cloudsc_fortran_atlas/expand_atlas_mod.F90
        ../cloudsc_fortran_atlas/validate_atlas_mod.F90
        ../cloudsc_gpu/cloudsc_gpu_scc_k_caching_mod.F90
        cloudsc_driver_scc_mod.F90
        dwarf_cloudsc_atlas.F90
    LIBS
        atlas_f
        pluto_f
        cloudsc-common-${prec}-atlas
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_K_CACHING
)

set_target_properties( dwarf-cloudsc-gpu-scc-k-caching-atlas-${prec}
        PROPERTIES
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module-dwarf-cloudsc-gpu-scc-k-caching-atlas"
)

ecbuild_add_test(
    TARGET dwarf-cloudsc-gpu-scc-k-caching-atlas-${prec}
    COMMAND bin/dwarf-cloudsc-gpu-scc-k-caching-atlas-${prec}
    ARGS 1 1000 128
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OMP 1
    ENVIRONMENT CLOUDSC_ATLAS_MULTIFIELD=0
)

ecbuild_add_test(
    TARGET dwarf-cloudsc-gpu-scc-k-caching-multifield-atlas-${prec}
    COMMAND bin/dwarf-cloudsc-gpu-scc-k-caching-atlas-${prec}
    ARGS 1 1000 128
    WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    OMP 1
    ENVIRONMENT CLOUDSC_ATLAS_MULTIFIELD=1
)
