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
    DESCRIPTION "Build optimized GPU version of CLOUDSC derived from CLAW" DEFAULT OFF
    CONDITION Serialbox_FOUND OR HDF5_FOUND
)

ecbuild_add_option( FEATURE CLOUDSC_GPU_SCC
    DESCRIPTION "Build optimized GPU version of CLOUDSC using SCC layout" DEFAULT OFF
    CONDITION Serialbox_FOUND OR HDF5_FOUND
)

ecbuild_add_option( FEATURE CLOUDSC_GPU_SCC_HOIST
    DESCRIPTION "Build optimized GPU version of CLOUDSC using SCC with hoisted temporary arrays" DEFAULT OFF
    CONDITION Serialbox_FOUND OR HDF5_FOUND
)


if( HAVE_CLOUDSC_GPU_CLAW )
    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-claw
        SOURCES
            dwarf_cloudsc_gpu.F90
            cloudsc_driver_gpu_claw_mod.F90
            cloudsc.claw.gpu.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_CLAW
    )
    if( HAVE_ACC AND TARGET OpenACC::OpenACC_Fortran )
        target_link_libraries( dwarf-cloudsc-gpu-claw PRIVATE OpenACC::OpenACC_Fortran )
        target_link_libraries( dwarf-cloudsc-gpu-claw PRIVATE cloudsc-common-lib-static )
    else()
        target_link_libraries( dwarf-cloudsc-gpu-claw PRIVATE cloudsc-common-lib )
    endif()

    if( HAVE_MPI AND TARGET MPI::MPI_Fortran )
        target_link_libraries( dwarf-cloudsc-gpu-claw PRIVATE MPI::MPI_Fortran )
    endif()
    if ( HAVE_HDF5 )
        target_include_directories( dwarf-cloudsc-gpu-claw PRIVATE ${HDF5_Fortran_INCLUDE_DIRS} )
    endif()
endif()


if( HAVE_CLOUDSC_GPU_SCC )
    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-scc
        SOURCES
            dwarf_cloudsc_gpu.F90
            cloudsc_driver_gpu_scc_mod.F90
            cloudsc_gpu_scc_mod.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC
    )
    if( HAVE_ACC AND TARGET OpenACC::OpenACC_Fortran )
        target_link_libraries( dwarf-cloudsc-gpu-scc PRIVATE OpenACC::OpenACC_Fortran )
        target_link_libraries( dwarf-cloudsc-gpu-scc PRIVATE cloudsc-common-lib-static )
    else()
        target_link_libraries( dwarf-cloudsc-gpu-scc PRIVATE cloudsc-common-lib )
    endif()

    if( HAVE_MPI AND TARGET MPI::MPI_Fortran )
        target_link_libraries( dwarf-cloudsc-gpu-scc PRIVATE MPI::MPI_Fortran )
    endif()
    if ( HAVE_HDF5 )
        target_include_directories( dwarf-cloudsc-gpu-scc PRIVATE ${HDF5_Fortran_INCLUDE_DIRS} )
    endif()

    # OpenMP5
    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-omp-scc
        SOURCES
            dwarf_cloudsc_gpu.F90
            cloudsc_driver_gpu_omp_scc_mod.F90
            cloudsc_gpu_omp_scc_mod.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC
    )
    if( HAVE_OMP AND TARGET OpenMP::OpenMP_Fortran )
        target_compile_definitions( dwarf-cloudsc-gpu-omp-scc PRIVATE CLOUDSC_GPU_OMP5)
        target_link_libraries( dwarf-cloudsc-gpu-omp-scc PRIVATE OpenACC::OpenACC_Fortran )
        target_link_libraries( dwarf-cloudsc-gpu-omp-scc PRIVATE cloudsc-common-lib-static )
    else()
        target_link_libraries( dwarf-cloudsc-gpu-omp-scc PRIVATE cloudsc-common-lib )
    endif()

    if( HAVE_MPI AND TARGET MPI::MPI_Fortran )
        target_link_libraries( dwarf-cloudsc-gpu-omp-scc PRIVATE MPI::MPI_Fortran )
    endif()
    if ( HAVE_HDF5 )
        target_include_directories( dwarf-cloudsc-gpu-omp-scc PRIVATE ${HDF5_Fortran_INCLUDE_DIRS} )
    endif()
endif()


if( HAVE_CLOUDSC_GPU_SCC_HOIST )
    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-scc-hoist
        SOURCES
            dwarf_cloudsc_gpu.F90
            cloudsc_driver_gpu_scc_hoist_mod.F90
            cloudsc_gpu_scc_hoist_mod.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_SCC_HOIST
    )
    if( HAVE_ACC AND TARGET OpenACC::OpenACC_Fortran )
        target_link_libraries( dwarf-cloudsc-gpu-scc-hoist PRIVATE OpenACC::OpenACC_Fortran )
        target_link_libraries( dwarf-cloudsc-gpu-scc-hoist PRIVATE cloudsc-common-lib-static )
    else()
        target_link_libraries( dwarf-cloudsc-gpu-scc-hoist PRIVATE cloudsc-common-lib )
    endif()

    if( HAVE_MPI AND TARGET MPI::MPI_Fortran )
        target_link_libraries( dwarf-cloudsc-gpu-scc-hoist PRIVATE MPI::MPI_Fortran )
    endif()
    if ( HAVE_HDF5 )
        target_include_directories( dwarf-cloudsc-gpu-scc-hoist PRIVATE ${HDF5_Fortran_INCLUDE_DIRS} )
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
