ecbuild_enable_fortran(REQUIRED MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module)

list(APPEND CLOUDSC_COMMON_SOURCES
  module/parkind1.F90
  module/yoecldp.F90
  module/yomcst.F90
  module/yoethf.F90
  module/yoephli.F90
  module/yomphyder.F90
  module/fcttre_mod.F90
  module/fccld_mod.F90
  module/routines.F90
  module/abor1.F90
  module/timer_mod.F90
  module/mycpu.c
  module/ec_pmon_mod.F90
  module/load_array_mod.F90
  module/validate_mod.F90
  module/cloudsc_global_state_mod.F90
  module/cloudsc_mpi_mod.F90
)


ecbuild_add_library( TARGET cloudsc-common-lib
  SOURCES ${CLOUDSC_COMMON_SOURCES}
  TYPE SHARED
  DEFINITIONS ${CLOUDSC_DEFINITIONS}
)
target_include_directories( cloudsc-common-lib PUBLIC $<INSTALL_INTERFACE:include> )
target_include_directories( cloudsc-common-lib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
target_include_directories( cloudsc-common-lib PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/module>)


if( HAVE_ACC )
  # ========================================================================
  # For OpenACC runs we can ensure global module variables from our utility
  # modules are copied to devices via `!$acc declare copyin(<module_vars>)`
  # However, this requries static linking as per this thread:
  # https://www.pgroup.com/userforum/viewtopic.php?t=7296
  # ========================================================================

  ecbuild_add_library( TARGET cloudsc-common-lib-static
    SOURCES ${CLOUDSC_COMMON_SOURCES}
    TYPE STATIC
    DEFINITIONS ${CLOUDSC_DEFINITIONS}
  )
  target_include_directories( cloudsc-common-lib-static PUBLIC $<INSTALL_INTERFACE:include> )
  target_include_directories( cloudsc-common-lib-static PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>)
  target_include_directories( cloudsc-common-lib-static PUBLIC $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/module>)
endif()

if( HAVE_OMP AND TARGET OpenMP::OpenMP_Fortran )
  target_link_libraries( cloudsc-common-lib PRIVATE OpenMP::OpenMP_Fortran )

  target_link_libraries( cloudsc-common-lib-static PRIVATE OpenMP::OpenMP_Fortran )
endif()

if( HAVE_MPI AND TARGET MPI::MPI_Fortran )
  target_sources( cloudsc-common-lib PRIVATE module/cloudsc_mpif.F90 )
  target_link_libraries( cloudsc-common-lib PRIVATE MPI::MPI_Fortran )

  target_sources( cloudsc-common-lib-static PRIVATE module/cloudsc_mpif.F90 )
  target_link_libraries( cloudsc-common-lib-static PRIVATE MPI::MPI_Fortran )
endif()

if ( HAVE_HDF5 )
  target_sources( cloudsc-common-lib PRIVATE module/hdf5_file_mod.F90 )
  target_include_directories( cloudsc-common-lib PRIVATE ${HDF5_Fortran_INCLUDE_DIRS} )
  target_link_libraries( cloudsc-common-lib PRIVATE ${HDF5_LIBRARIES} )

  target_sources( cloudsc-common-lib-static PRIVATE module/hdf5_file_mod.F90 )
  target_include_directories( cloudsc-common-lib-static PRIVATE ${HDF5_Fortran_INCLUDE_DIRS} )
  target_link_libraries( cloudsc-common-lib-static PRIVATE ${HDF5_LIBRARIES} )
endif()

if( HAVE_ACC AND TARGET OpenACC::OpenACC_Fortran )
  target_link_libraries( cloudsc-common-lib PRIVATE OpenACC::OpenACC_Fortran )

  target_link_libraries( cloudsc-common-lib-static PRIVATE OpenACC::OpenACC_Fortran )
endif()

if( TARGET Serialbox::Serialbox_Fortran )
  target_link_libraries( cloudsc-common-lib PRIVATE Serialbox::Serialbox_Fortran )

  target_link_libraries( cloudsc-common-lib-static PRIVATE Serialbox::Serialbox_Fortran )
endif()
