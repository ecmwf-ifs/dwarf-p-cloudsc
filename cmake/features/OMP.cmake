if( HAVE_OMP )

  if( CMAKE_Fortran_COMPILER_ID MATCHES "Cray" )
    # Workaround for Linker issue with Cray compiler, see
    # https://gitlab.kitware.com/cmake/cmake/-/issues/24402
    set( OMP_LINK_OPTIONS "-fopenmp" )
  else()
    set( OMP_LINK_OPTIONS )
  endif()


  if( NOT DEFINED HAVE_OMP_TARGET_TEAMS_DISTRIBUTE )

    try_compile(
      HAVE_OMP_TARGET_TEAMS_DISTRIBUTE
      ${CMAKE_CURRENT_BINARY_DIR}
      ${PROJECT_SOURCE_DIR}/cmake/features/OMP/test_omp_target_teams_distribute.F90
      LINK_LIBRARIES OpenMP::OpenMP_Fortran
      LINK_OPTIONS ${OMP_LINK_OPTIONS}
      OUTPUT_VARIABLE _HAVE_OMP_TARGET_TEAMS_DISTRIBUTE_OUTPUT
    )

    ecbuild_debug_var( HAVE_OMP_TARGET_TEAMS_DISTRIBUTE )
    ecbuild_debug_var( _HAVE_OMP_TARGET_TEAMS_DISTRIBUTE_OUTPUT )

  endif()

  if( NOT DEFINED HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL )

    try_compile(
      HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL
      ${CMAKE_CURRENT_BINARY_DIR}
      ${PROJECT_SOURCE_DIR}/cmake/features/OMP/test_omp_target_loop_construct_bind_parallel.F90
      LINK_LIBRARIES OpenMP::OpenMP_Fortran
      LINK_OPTIONS ${OMP_LINK_OPTIONS}
      OUTPUT_VARIABLE _HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL_OUTPUT
    )

    ecbuild_debug_var( HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL )
    ecbuild_debug_var( _HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL_OUTPUT )

  endif()

  if( NOT DEFINED HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD )

    try_compile(
      HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD
      ${CMAKE_CURRENT_BINARY_DIR}
      ${PROJECT_SOURCE_DIR}/cmake/features/OMP/test_omp_target_loop_construct_bind_thread.F90
      LINK_LIBRARIES OpenMP::OpenMP_Fortran
      LINK_OPTIONS ${OMP_LINK_OPTIONS}
      OUTPUT_VARIABLE _HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD_OUTPUT
    )

    ecbuild_debug_var( HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD )
    ecbuild_debug_var( _HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD_OUTPUT )

  endif()

  if( HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL OR HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD )
    set( HAVE_OMP_TARGET_LOOP_CONSTRUCT ON CACHE BOOL "Support for OpenMP target teams loop" )
  else()
    set( HAVE_OMP_TARGET_LOOP_CONSTRUCT OFF CACHE BOOL "Support for OpenMP target teams loop" )
  endif()

  if( HAVE_OMP_TARGET_TEAMS_DISTRIBUTE OR HAVE_OMP_TARGET_LOOP_CONSTRUCT )
    set( HAVE_OMP_TARGET ON CACHE BOOL "Support for OpenMP target" )
  else()
    set( HAVE_OMP_TARGET OFF CACHE BOOL "Support for OpenMP target" )
  endif()

  ecbuild_debug_var( HAVE_OMP_TARGET_LOOP_CONSTRUCT )
  ecbuild_debug_var( HAVE_OMP_TARGET )

endif()
