# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# define minimum version of cmake required
cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

find_package( ecbuild REQUIRED )

# define the project
project( dwarf-p-cloudsc LANGUAGES C Fortran )

include( cmake/compat.cmake )
if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  ecbuild_add_fortran_flags("-ffree-line-length-none")
  if( CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10.0")
    ecbuild_add_fortran_flags("-fallow-argument-mismatch")
  endif()
endif()

### OpenACC
if( NOT DEFINED ENABLE_ACC OR ENABLE_ACC )
    find_package( OpenACC COMPONENTS Fortran C )
    # Ensure targets are available
    if( OpenACC_Fortran_FOUND AND NOT TARGET OpenACC::OpenACC_Fortran )
        ecbuild_critical( "OpenACC::OpenACC_Fortran should have been defined as a CMake target" )
    endif()
    if( OpenACC_C_FOUND AND NOT TARGET OpenACC::OpenACC_C )
        ecbuild_critical( "OpenACC::OpenACC_C should have been defined as a CMake target" )
    endif()
endif()
ecbuild_add_option( FEATURE ACC
    DESCRIPTION "OpenACC" DEFAULT ON
    CONDITION OpenACC_C_FOUND OR OpenACC_Fortran_FOUND )

### OpenMP
if( NOT DEFINED ENABLE_OMP OR ENABLE_OMP )
    find_package( OpenMP COMPONENTS Fortran C )
    # Ensure targets are available
    if( OpenMP_Fortran_FOUND AND NOT TARGET OpenMP::OpenMP_Fortran )
        ecbuild_critical( "OpenMP::OpenMP_Fortran should have been defined as a CMake target" )
    endif()
    if( OpenMP_C_FOUND AND NOT TARGET OpenMP::OpenMP_C )
        ecbuild_critical( "OpenMP::OpenMP_C should have been defined as a CMake target" )
    endif()
endif()
ecbuild_add_option( FEATURE OMP
    DESCRIPTION "OpenMP" DEFAULT ON
    CONDITION OpenMP_C_FOUND AND OpenMP_Fortran_FOUND )

### MPI
if( NOT DEFINED ENABLE_MPI OR ENABLE_MPI )
    ecbuild_find_mpi( COMPONENTS Fortran )
endif()
if( MPI_Fortran_FOUND AND ENABLE_MPI )
    # Ensure targets are available
    if( MPI_Fortran_FOUND AND NOT TARGET MPI::MPI_Fortran )
        ecbuild_critical( "MPI::MPI_Fortran should have been defined as a CMake target" )
    endif()
    ecbuild_add_option( FEATURE MPI
        DESCRIPTION "MPI" DEFAULT OFF
        CONDITION MPI_Fortran_FOUND )
    list(APPEND CLOUDSC_DEFINITIONS HAVE_MPI ${MPI_DEFINITIONS})
endif()

### HDF5
ecbuild_add_option( FEATURE HDF5
    DESCRIPTION "Use HDF5 to read input and reference data"
    REQUIRED_PACKAGES "HDF5 COMPONENTS Fortran"
    CONDITION NOT DEFINED ENABLE_HDF5 OR ENABLE_HDF5
    DEFAULT ON )
if( HDF5_FOUND )
    list(APPEND CLOUDSC_DEFINITIONS HAVE_HDF5 ${HDF5_Fortran_DEFINITIONS})
endif()

# Add Serialbox utility package for platform-agnostic file I/O
ecbuild_add_option( FEATURE SERIALBOX
    DESCRIPTION "Use Serialbox to read input and reference data"
    REQUIRED_PACKAGES "Serialbox"
    CONDITION (NOT DEFINED ENABLE_SERIALBOX OR ENABLE_SERIALBOX) AND NOT HDF5_FOUND
    DEFAULT OFF )
if( SERIALBOX_FOUND )
    list(APPEND CLOUDSC_DEFINITIONS HAVE_SERIALBOX)
endif()

ecbuild_find_package( NAME loki )

# Add option for single-precision builds
ecbuild_add_option( FEATURE SINGLE_PRECISION
    DESCRIPTION "Build CLOUDSC in single precision" DEFAULT OFF
)
if( HAVE_SINGLE_PRECISION )
    list(APPEND CLOUDSC_DEFINITIONS SINGLE)
endif()

# build executables
add_subdirectory(src)

# documentation
find_package( Latex )
ecbuild_add_option(FEATURE DOCS
                   DESCRIPTION "Documentation"
                   DEFAULT OFF
                   CONDITION Latex_FOUND)


if(HAVE_DOCS)
  add_subdirectory(doc)
endif()

# finalize
ecbuild_install_project(NAME dwarf-p-cloudsc)

# print summary
ecbuild_print_summary()
