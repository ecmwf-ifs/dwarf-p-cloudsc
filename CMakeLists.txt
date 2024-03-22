# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# define minimum version of cmake required
cmake_minimum_required( VERSION 3.17 FATAL_ERROR )

find_package( ecbuild REQUIRED )

# define the project
project( dwarf-p-cloudsc LANGUAGES C Fortran )

include( cmake/compat.cmake )
if( CMAKE_Fortran_COMPILER_ID MATCHES "GNU")
  ecbuild_add_fortran_flags("-ffree-line-length-none")
  if( CMAKE_Fortran_COMPILER_VERSION VERSION_GREATER_EQUAL "10.0")
    ecbuild_add_fortran_flags("-fallow-argument-mismatch")
  endif()
  if( CMAKE_Fortran_COMPILER_VERSION VERSION_LESS "11.0")
    set( ENABLE_ACC OFF )
  endif()
endif()

include( cloudsc_compile_options )

### OpenACC
if( ${CMAKE_VERSION} VERSION_LESS "3.25" AND (NOT DEFINED ENABLE_ACC OR ENABLE_ACC ) )
    # Incredibly inconvenient: FindOpenACC does _not_ set OpenACC_FOUND, only
    # the language-specific components OpenACC_Fortran_FOUND and OpenACC_C_FOUND.
    # This means, even internally CMake considers OpenACC as not found.
    # (See eg get_property(... GLOBAL PROPERTY PACKAGES_NOT_FOUND))
    # Therefore, we search for OpenACC, set OpenACC_FOUND ourselves according to
    # the result, and then, trigger a second find_package via ecbuild_add_option.
    # This then conveniently takes the previously set OpenACC_FOUND into account
    # and rectifies CMake's internal bookkeeping in the process.
    # This has been fixed in CMake 3.25
    find_package( OpenACC )
    if( OpenACC_Fortran_FOUND AND OpenACC_C_FOUND )
        set( OpenACC_FOUND ON )
    endif()
endif()
ecbuild_add_option( FEATURE ACC
    DESCRIPTION "OpenACC" DEFAULT ON
    REQUIRED_PACKAGES "OpenACC" )

include(CheckLanguage)
check_language( CUDA )
ecbuild_add_option( FEATURE CUDA
    DESCRIPTION "CUDA" DEFAULT OFF
    CONDITION CMAKE_CUDA_COMPILER )
if( HAVE_CUDA )
    enable_language( CUDA )
endif()

ecbuild_add_option( FEATURE HIP
	DESCRIPTION "HIP" DEFAULT OFF
	REQUIRED_PACKAGES "hip" )
if ( HAVE_HIP )
  if(NOT DEFINED ROCM_PATH)
    if(DEFINED ENV{ROCM_PATH})
      set(ROCM_PATH $ENV{ROCM_PATH} CACHE PATH "Path to which ROCM has been installed")
    else()
      set(ROCM_PATH "/opt/rocm" CACHE PATH "Path to which ROCM has been installed")
    endif()
  endif()
  find_package(hip REQUIRED)
endif()

ecbuild_add_option( FEATURE SYCL
	DESCRIPTION "SYCL" DEFAULT OFF
	REQUIRED_PACKAGES "IntelSYCL")

### OpenMP
ecbuild_add_option( FEATURE OMP
    DESCRIPTION "OpenMP" DEFAULT ON
    REQUIRED_PACKAGES "OpenMP COMPONENTS Fortran C" )
include(features/OMP)

### MPI
ecbuild_add_option( FEATURE MPI
    DESCRIPTION "MPI" DEFAULT OFF
    REQUIRED_PACKAGES "MPI COMPONENTS Fortran" )
if( HAVE_MPI )
    list(APPEND CLOUDSC_DEFINITIONS HAVE_MPI )
endif()

### HDF5
ecbuild_add_option( FEATURE HDF5
    DESCRIPTION "Use HDF5 to read input and reference data"
    REQUIRED_PACKAGES "HDF5 COMPONENTS Fortran C"
    DEFAULT ON )
if( HAVE_HDF5 )
    list(APPEND CLOUDSC_DEFINITIONS HAVE_HDF5 )
endif()

# Add Serialbox utility package for platform-agnostic file I/O
ecbuild_add_option( FEATURE SERIALBOX
    DESCRIPTION "Use Serialbox to read input and reference data"
    REQUIRED_PACKAGES "Serialbox"
    CONDITION NOT HAVE_HDF5
    DEFAULT OFF )
if( HAVE_SERIALBOX )
    list(APPEND CLOUDSC_DEFINITIONS HAVE_SERIALBOX)
endif()

# Add field_api library to manage data fields
ecbuild_add_option( FEATURE FIELD_API
    DESCRIPTION "Use field_api to manage GPU data offload and copyback"
    REQUIRED_PACKAGES "field_api"
    CONDITION HAVE_CUDA
    DEFAULT ON )

ecbuild_find_package( NAME loki )
ecbuild_find_package( NAME atlas )

# Add option for single-precision builds
ecbuild_add_option( FEATURE SINGLE_PRECISION
    DESCRIPTION "Build CLOUDSC in single precision" DEFAULT OFF
)
set(prec dp)
if( HAVE_SINGLE_PRECISION )
    list(APPEND CLOUDSC_DEFINITIONS SINGLE)
    set(prec sp)
endif()

# build executables
add_subdirectory(src)

# documentation
ecbuild_add_option(FEATURE DOCS
                   DESCRIPTION "Documentation"
                   REQUIRED_PACKAGES "Latex"
                   DEFAULT OFF)


if(HAVE_DOCS)
  add_subdirectory(doc)
endif()

# finalize
ecbuild_install_project(NAME dwarf-p-cloudsc)

# print summary
ecbuild_print_summary()
