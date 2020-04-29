# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_GPU_NAIVE
    DESCRIPTION "Build the naive GPU version CLOUDSC using `!$acc kernels`" DEFAULT ON
)

if( HAVE_CLOUDSC_GPU_NAIVE )
    ecbuild_add_executable(
        TARGET dwarf-cloudsc-gpu-naive
        SOURCES
            dwarf_cloudsc_gpu_naive.F90
            cloudsc_global_state_mod.F90
            cloudsc_driver_gpu_naive_mod.F90
            cloudsc_acc_kernels_mod.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    target_link_libraries( dwarf-cloudsc-gpu-naive PRIVATE cloudsc-common-lib )

    if( TARGET OpenACC::OpenACC_Fortran )
        target_link_libraries( dwarf-cloudsc-gpu-naive PRIVATE OpenACC::OpenACC_Fortran )
    endif()

    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )
endif()
