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

if( HAVE_CLOUDSC_GPU_KERNELS )
    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-kernels
        SOURCES
            dwarf_cloudsc_gpu.F90
            cloudsc_driver_gpu_kernels_mod.F90
            cloudsc_acc_kernels_mod.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_KERNELS
    )
    target_link_libraries( dwarf-cloudsc-gpu-kernels PRIVATE cloudsc-common-lib )

    if( HAVE_ACC AND TARGET OpenACC::OpenACC_Fortran )
        target_link_libraries( dwarf-cloudsc-gpu-kernels PRIVATE OpenACC::OpenACC_Fortran )
    endif()

    if( HAVE_MPI AND TARGET MPI::MPI_Fortran )
        target_link_libraries( dwarf-cloudsc-gpu-kernels PRIVATE MPI::MPI_Fortran )
    endif()
endif()


if( HAVE_CLOUDSC_GPU_CLAW )
    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-claw
        SOURCES
            dwarf_cloudsc_gpu.F90
            cloudsc_driver_gpu_claw_mod.F90
            cloudsc.claw.gpu.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS} CLOUDSC_GPU_CLAW
    )
    target_link_libraries( dwarf-cloudsc-gpu-claw PRIVATE cloudsc-common-lib )

    if( HAVE_ACC AND TARGET OpenACC::OpenACC_Fortran )
        target_link_libraries( dwarf-cloudsc-gpu-claw PRIVATE OpenACC::OpenACC_Fortran )
    endif()

    if( HAVE_MPI AND TARGET MPI::MPI_Fortran )
        target_link_libraries( dwarf-cloudsc-gpu-claw PRIVATE MPI::MPI_Fortran )
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
