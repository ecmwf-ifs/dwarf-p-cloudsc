###############################################################################
# (C) Copyright 1996-2015 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.
###############################################################################

# define minimum version of cmake required
cmake_minimum_required( VERSION 3.12 FATAL_ERROR )

set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild/cmake" )
set( CMAKE_BUILD_TYPE None )

set(ENABLE_OS_TYPES_TEST       OFF   CACHE BOOL " ")
set(ENABLE_OS_FUNCTIONS_TEST   OFF   CACHE BOOL " ")
set(ENABLE_LARGE_FILE_SUPPORT  OFF   CACHE BOOL " ")
set(ENABLE_OS_ENDINESS_TEST    OFF   CACHE BOOL " ")

# define the project
project(dwarf-P-cloudMicrophysics-IFSScheme VERSION 1.1 LANGUAGES C Fortran )

# include ecbuild required commands
include(ecbuild_system NO_POLICY_SCOPE)
ecbuild_requires_macro_version(2.10.1)
ecbuild_declare_project()

include(CompilerFlags)

find_package( OpenMP COMPONENTS Fortran C )

# Add Serialbox utility package for platform-agnostic file I/O
find_package( Serialbox )
if( Serialbox_FOUND )
    list(APPEND CLOUDSC_DEFINITIONS HAVE_SERIALBOX)
endif()

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
ecbuild_install_project(NAME dwarf-P-cloudMicrophysics-IFSScheme)

# print summary
ecbuild_print_summary()
