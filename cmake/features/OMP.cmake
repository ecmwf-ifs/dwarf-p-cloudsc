if( HAVE_OMP )

    try_compile(
        HAVE_OMP_TARGET_TEAMS_DISTRIBUTE
        ${CMAKE_CURRENT_BINARY_DIR}
        ${PROJECT_SOURCE_DIR}/cmake/features/OMP/test_omp_target_teams_distribute.F90
        LINK_LIBRARIES OpenMP::OpenMP_Fortran
        OUTPUT_VARIABLE _HAVE_OMP_TARGET_TEAMS_DISTRIBUTE_OUTPUT
    )

    ecbuild_debug_var( HAVE_OMP_TARGET_TEAMS_DISTRIBUTE )
    ecbuild_debug_var( _HAVE_OMP_TARGET_TEAMS_DISTRIBUTE_OUTPUT )

    try_compile(
        HAVE_OMP_TARGET_LOOP_CONSTRUCT
        ${CMAKE_CURRENT_BINARY_DIR}
        ${PROJECT_SOURCE_DIR}/cmake/features/OMP/test_omp_target_loop_construct.F90
        LINK_LIBRARIES OpenMP::OpenMP_Fortran
        OUTPUT_VARIABLE _HAVE_OMP_TARGET_LOOP_CONSTRUCT_OUTPUT
    )

    ecbuild_debug_var( HAVE_OMP_TARGET_LOOP_CONSTRUCT )
    ecbuild_debug_var( _HAVE_OMP_TARGET_LOOP_CONSTRUCT_OUTPUT )

endif()
