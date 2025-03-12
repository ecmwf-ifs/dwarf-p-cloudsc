# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# define minimum version of cmake required
# 3.21 is required for HIP support
# 3.23.5 is required for SYCL support
cmake_minimum_required( VERSION 3.23.5 FATAL_ERROR )

find_package( ecbuild REQUIRED )

# define the project
project( dwarf-p-cloudsc LANGUAGES C CXX Fortran )

ecbuild_enable_fortran( REQUIRED NO_MODULE_DIRECTORY )

include( cloudsc_compile_options )

### Precision selection

set( CLOUDSC_DEFINITIONS )

ecbuild_add_option( FEATURE DOUBLE_PRECISION
                    DEFAULT ON
                    DESCRIPTION "Build double precision CLOUDSC binaries" )

ecbuild_add_option( FEATURE SINGLE_PRECISION
                    DEFAULT ON
                    DESCRIPTION "Build single precision CLOUDSC binaries" )

### Programming models

ecbuild_add_option( FEATURE OMP
                    DESCRIPTION "Support for thread-level parallelism or GPUs using OpenMP"
                    DEFAULT ON
                    REQUIRED_PACKAGES "OpenMP COMPONENTS Fortran C CXX" )
include(features/OMP)

ecbuild_add_option( FEATURE MPI
                    DESCRIPTION "Support for task-level parallelism using MPI"
                    DEFAULT OFF
                    REQUIRED_PACKAGES "MPI COMPONENTS Fortran" )
if( HAVE_MPI )
  list( APPEND CLOUDSC_DEFINITIONS HAVE_MPI )
endif()

include(features/ACC)
ecbuild_add_option( FEATURE ACC
                    DEFAULT ON
                    DESCRIPTION "Support for using GPUs with OpenACC"
                    REQUIRED_PACKAGES "OpenACC COMPONENTS Fortran" )

include(CheckLanguage)
check_language( CUDA )
ecbuild_add_option( FEATURE CUDA
                    DESCRIPTION "Support for using GPUs with CUDA"
                    DEFAULT ON
                    CONDITION CMAKE_CUDA_COMPILER
                    REQUIRED_PACKAGES "CUDAToolkit" )
include(features/CUDA)

check_language( HIP )
ecbuild_add_option( FEATURE HIP
	                  DESCRIPTION "Support for using GPUs with HIP"
                    DEFAULT ON
                    CONDITION CMAKE_HIP_COMPILER
                    REQUIRED_PACKAGES "hip" )
if( HAVE_HIP )
  enable_language( HIP )
endif()

ecbuild_add_option( FEATURE SYCL
	                  DESCRIPTION "Support for using GPUs with SYCL"
                    DEFAULT ON
                  	REQUIRED_PACKAGES "IntelSYCL" )

### Link options

ecbuild_add_option( FEATURE GPU_STATIC
                    DEFAULT ON
                    DESCRIPTION "Compile GPU targets using static libraries" )

### I/O libraries

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

### Enable/disable specific variants

# Add loki source-to-source translation toolchain to generate code variants
ecbuild_add_option( FEATURE LOKI
                    DESCRIPTION "Support for generating GPU variants with Loki"
                    DEFAULT ON
                    REQUIRED_PACKAGES "loki" )

# Add field_api library to manage data fields
ecbuild_add_option( FEATURE FIELD_API
                    DESCRIPTION "Use field_api to manage GPU data offload and copyback"
                    DEFAULT ON
                    REQUIRED_PACKAGES "field_api" )

ecbuild_add_option( FEATURE FIELD_API_DISABLE_MAPPED_MEMORY
                    DESCRIPTION "Disable the use of ACC mapped memory in Field API objects"
                    CONDITION HAVE_FIELD_API
                    DEFAULT OFF )
if( HAVE_FIELD_API_DISABLE_MAPPED_MEMORY )
  list(APPEND CLOUDSC_DEFINITIONS FIELD_API_DISABLE_MAPPED_MEMORY)
endif()

# Add atlas library to manage data fields
ecbuild_add_option( FEATURE ATLAS
                    DESCRIPTION "Use Atlas to manage GPU data offload and copyback"
                    DEFAULT ON
                    REQUIRED_PACKAGES "atlas" )

# A variant driven by Python
# NB: This requires Python3 to install f90wrap and its dependencies, and has only been
#     tested successfully with GNU and Intel compilers
find_package( Python3 COMPONENTS Interpreter QUIET )
ecbuild_add_option( FEATURE PYIFACE
                    DESCRIPTION "Build a (slightly customized) cloudsc-fortran variant with Python interface"
                    DEFAULT OFF
                    CONDITION HAVE_HDF5 AND Python3_FOUND AND CMAKE_Fortran_COMPILER_ID MATCHES "GNU|Intel" )

ecbuild_add_option( FEATURE PYIFACE_BINARY
                    DESCRIPTION "Build the debug-oriented binary for Python interface of (customized) cloudsc-fortran"
                    DEFAULT OFF
                    CONDITION HAVE_PYIFACE )

# A pure Python variant, created from Loki transpilation
ecbuild_add_option( FEATURE PYTHON_F2PY
                    DESCRIPTION "Build the pure Python variant from Loki transpilation"
                    DEFAULT OFF
                    CONDITION HAVE_HDF5 AND Python3_FOUND AND CMAKE_Fortran_COMPILER_ID MATCHES "GNU|Intel" )

### Input data

if( HAVE_HDF5 )
  file( CREATE_LINK
        "${PROJECT_SOURCE_DIR}/config-files/input.h5"
        "${CMAKE_BINARY_DIR}/input.h5" )
  file( CREATE_LINK
        "${PROJECT_SOURCE_DIR}/config-files/reference.h5"
        "${CMAKE_BINARY_DIR}/reference.h5" )
elseif( HAVE_SERIALBOX )
  file( CREATE_LINK
        "${PROJECT_SOURCE_DIR}/data"
        "${CMAKE_BINARY_DIR}/data" )
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
