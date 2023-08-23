# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

ecbuild_enable_fortran(REQUIRED MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/module)

list(APPEND CLOUDSC_COMMON_SOURCES
    module/parkind1.F90
    module/yoecldp.F90
    module/yomcst.F90
    module/yoethf.F90
    module/yoephli.F90
    module/yomphyder.F90
    module/routines.F90
    module/abor1.F90
    module/timer_mod.F90
    module/mycpu.c
    module/ec_pmon_mod.F90
    module/oml_mod.F90
    module/file_io_mod.F90
    module/expand_mod.F90
    module/validate_mod.F90
    module/cloudsc_global_state_mod.F90
    module/cloudsc_mpi_mod.F90
)

if(NOT ENABLE_USE_STMT_FUNC)
    list(APPEND CLOUDSC_COMMON_SOURCES
        module/fcttre_mod.F90
        module/fccld_mod.F90
    )
endif()

list(APPEND CLOUDSC_CUDA_SOURCES
    module/fccld_mod.cuf.F90
    module/fcttre_mod.cuf.F90
    module/yoethf.cuf.F90
    module/yomcst.cuf.F90
)

if( HAVE_CUDA )
    # ========================================================================
    # Compile CUDA fortran files with -cuda.
    #
    # This is necessary since CMake's CUDA languages does not natively
    # understand CUDA-Fortran (.cuf) yet. So we simply emulate .cuf with
    # .cuf.F90 files and force the compile flags here.
    # ========================================================================
    cloudsc_add_compile_options(
        SOURCES ${CLOUDSC_CUDA_SOURCES}
        FLAGS "-cuda -gpu=maxregcount:128"
    )

    # Add CUDA-specific flags to the library if enabled
    list(APPEND CLOUDSC_COMMON_SOURCES ${CLOUDSC_CUDA_SOURCES} )
endif()


if( HAVE_ACC )
    # ========================================================================
    # For OpenACC runs we can ensure global module variables from our utility
    # modules are copied to devices via `!$acc declare copyin(<module_vars>)`
    # However, this requries static linking as per this thread:
    # https://www.pgroup.com/userforum/viewtopic.php?t=7296
    # ========================================================================
    set( LIBRARY_TYPE STATIC )
else()
    set( LIBRARY_TYPE SHARED )
endif()

ecbuild_add_library( TARGET cloudsc-common-lib
    TYPE ${LIBRARY_TYPE}
    DEFINITIONS
        ${CLOUDSC_DEFINITIONS}
	    $<${HAVE_FIELD_API}:USE_FIELD_API>
    SOURCES
        ${CLOUDSC_COMMON_SOURCES}
        $<${HAVE_MPI}:module/cloudsc_mpif.F90>
        $<${HAVE_HDF5}:module/hdf5_file_mod.F90>
        $<${HAVE_FIELD_API}:module/cloudsc_field_state_mod.F90>
    PRIVATE_INCLUDES
        $<${HAVE_HDF5}:${HDF5_Fortran_INCLUDE_DIRS}>
    PUBLIC_INCLUDES
        $<INSTALL_INTERFACE:include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
        $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/module>
    PUBLIC_LIBS
        $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
        $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
        $<${HAVE_MPI}:MPI::MPI_Fortran>
        $<${HAVE_HDF5}:hdf5::hdf5_fortran>
        $<${HAVE_SERIALBOX}:Serialbox::Serialbox_Fortran>
        $<${HAVE_FIELD_API}:field_api>
)

if( HAVE_CUDA )
    target_link_options( cloudsc-common-lib PUBLIC "-cuda" )
endif()
