# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

### Source files for common lib

set(CLOUDSC_COMMON_SOURCES
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

set(CLOUDSC_CUDA_SOURCES
    module/fccld_mod.cuf.F90
    module/fcttre_mod.cuf.F90
    module/yoethf.cuf.F90
    module/yomcst.cuf.F90
)

### Helper macro to build common lib variant

macro( cloudsc_add_common_lib )

    set( options "" )
    set( oneValueArgs SUFFIX TYPE )
    set( multiValueArgs DEFINITIONS SOURCES LIBS )

    cmake_parse_arguments( COMMON_LIB "${options}" "${oneValueArgs}" "${multiValueArgs}" ${ARGN} )
    if( COMMON_LIB_UNPARSED_ARGUMENTS )
        ecbuild_critical( "Unknown arguments given to add_cloudsc_common_lib(): \"${COMMON_LIB_UNPARSED_ARGUMENTS}\"")
    endif()

    ecbuild_add_library(
        TARGET cloudsc-common-${prec}-${COMMON_LIB_SUFFIX}
        TYPE ${COMMON_LIB_TYPE}
        PUBLIC_DEFINITIONS
            ${CLOUDSC_DEFINITIONS}
            $<$<STREQUAL:"${prec}","sp">:SINGLE>
            ${COMMON_LIB_DEFINITIONS}
        SOURCES
            ${CLOUDSC_COMMON_SOURCES}
            ${COMMON_LIB_SOURCES}
            $<${HAVE_HDF5}:module/hdf5_file_mod.F90>
            $<${HAVE_MPI}:module/cloudsc_mpif.F90>
        PUBLIC_LIBS
            $<${HAVE_HDF5}:hdf5::hdf5_fortran>
            $<${HAVE_SERIALBOX}:Serialbox::Serialbox_Fortran>
            $<${HAVE_MPI}:MPI::MPI_Fortran>
            ${COMMON_LIB_LIBS}
        PUBLIC_INCLUDES
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_BINARY_DIR}/module-${COMMON_LIB_SUFFIX}>
    )

    set_target_properties( cloudsc-common-${prec}-${COMMON_LIB_SUFFIX}
        PROPERTIES
            Fortran_MODULE_DIRECTORY "${CMAKE_CURRENT_BINARY_DIR}/module-${COMMON_LIB_SUFFIX}"
    )

    if( CMAKE_C_COMPILER_ID MATCHES "Cray" AND HAVE_OMP )
        if( "${COMMON_LIB_LIBS}" MATCHES "OpenMP::OpenMP_Fortran" )
            # Workaround for Linker issue with Cray compiler, see
            # https://gitlab.kitware.com/cmake/cmake/-/issues/24402
            target_link_options( cloudsc-common-${prec}-${COMMON_LIB_SUFFIX} PUBLIC "-fopenmp" )
        endif()
    endif()

endmacro()

### Link type selection

if( BUILD_SHARED_LIBS )
    set( LINK_TYPE SHARED )
else()
    set( LINK_TYPE STATIC )
endif()

if( HAVE_GPU_STATIC )
    set( GPU_LINK_TYPE STATIC )
else()
    set( GPU_LINK_TYPE ${LINK_TYPE} )
endif()

### Common lib with OpenMP

cloudsc_add_common_lib(
    SUFFIX lib
    TYPE ${LINK_TYPE}
    LIBS
        $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
)

### Common lib for field_api-managed data structures

if( HAVE_FIELD_API )
    cloudsc_add_common_lib(
        SUFFIX field_api
        TYPE ${GPU_LINK_TYPE}
        DEFINITIONS USE_FIELD_API
        SOURCES
            module/cloudsc_field_state_mod.F90
            module/cloudsc_state_type_mod.F90
            module/cloudsc_aux_type_mod.F90
            module/cloudsc_flux_type_mod.F90
        LIBS
            $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
            $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
            field_api_${prec}
    )
endif()

### Common lib for Atlas-managed data structures

if( HAVE_ATLAS )
    # NB: This does not contain the global state and expand
    #     modules at the moment
    cloudsc_add_common_lib(
        SUFFIX atlas
        TYPE ${GPU_LINK_TYPE}
        LIBS
            $<${HAVE_OMP}:OpenMP::OpenMP_Fortran>
            $<${HAVE_ACC}:OpenACC::OpenACC_Fortran>
    )
endif()

### Common lib for OpenMP-target offload

if( HAVE_OMP AND HAVE_OMP_TARGET )

    cloudsc_add_common_lib(
        SUFFIX omp
        TYPE ${GPU_LINK_TYPE}
        LIBS
            OpenMP::OpenMP_Fortran
    )

    if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
        target_compile_options( cloudsc-common-${prec}-omp PUBLIC "-mp=gpu" )
    endif()

endif()

### Common lib for OpenACC support

if( HAVE_ACC )

    cloudsc_add_common_lib(
        SUFFIX acc
        TYPE ${GPU_LINK_TYPE}
        LIBS
            OpenACC::OpenACC_Fortran
    )

endif()

### Common lib for CUDA targets

if( HAVE_CUDA )
    # ========================================================================
    # Compile CUDA fortran files with -cuda.
    #
    # This is necessary since CMake's CUDA language does not natively
    # understand CUDA-Fortran (.cuf) yet. So we simply emulate .cuf with
    # .cuf.F90 files and force the compile flags here.
    # ========================================================================
    if( CMAKE_Fortran_COMPILER_ID MATCHES "NVHPC" )
        cloudsc_add_compile_options(
            SOURCES ${CLOUDSC_CUDA_SOURCES}
            FLAGS "-cuda -gpu=maxregcount:128"
        )
    endif()

    cloudsc_add_common_lib(
        SUFFIX cuda
        TYPE ${GPU_LINK_TYPE}
        SOURCES
            ${CLOUDSC_CUDA_SOURCES}
        LIBS
            OpenACC::OpenACC_Fortran
            CUDA::cudart
    )

    target_link_options( cloudsc-common-${prec}-cuda
        PUBLIC "-cuda" )

endif()
