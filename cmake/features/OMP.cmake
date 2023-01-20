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
        HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL
        ${CMAKE_CURRENT_BINARY_DIR}
        ${PROJECT_SOURCE_DIR}/cmake/features/OMP/test_omp_target_loop_construct_bind_parallel.F90
        LINK_LIBRARIES OpenMP::OpenMP_Fortran
        OUTPUT_VARIABLE _HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL_OUTPUT
    )

    ecbuild_debug_var( HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL )
    ecbuild_debug_var( _HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL_OUTPUT )

    try_compile(
        HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD
        ${CMAKE_CURRENT_BINARY_DIR}
        ${PROJECT_SOURCE_DIR}/cmake/features/OMP/test_omp_target_loop_construct_bind_thread.F90
        LINK_LIBRARIES OpenMP::OpenMP_Fortran
        OUTPUT_VARIABLE _HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD_OUTPUT
    )

    ecbuild_debug_var( HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD )
    ecbuild_debug_var( _HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD_OUTPUT )

    if( HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL OR HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD )
        set( HAVE_OMP_TARGET_LOOP_CONSTRUCT ON CACHE BOOL "OpenMP target teams loop is supported" )
    else()
        set( HAVE_OMP_TARGET_LOOP_CONSTRUCT OFF CACHE BOOL "OpenMP target teams loop is not supported" )
    endif()

    ecbuild_debug_var( HAVE_OMP_TARGET_LOOP_CONSTRUCT )

endif()
