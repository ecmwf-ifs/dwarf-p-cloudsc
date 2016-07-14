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
cmake_minimum_required(VERSION 2.8.4 FATAL_ERROR)

set( CMAKE_MODULE_PATH "${CMAKE_CURRENT_SOURCE_DIR}/cmake" ${CMAKE_MODULE_PATH} "${CMAKE_CURRENT_SOURCE_DIR}/../ecbuild/cmake" )
set( CMAKE_BUILD_TYPE None )

# define the project
project(dwarf3 )

set(ENABLE_OS_TYPES_TEST       OFF   CACHE BOOL " ")
set(ENABLE_OS_FUNCTIONS_TEST   OFF   CACHE BOOL " ")
set(ENABLE_LARGE_FILE_SUPPORT  OFF   CACHE BOOL " ")
set(ENABLE_OS_ENDINESS_TEST    OFF   CACHE BOOL " ")

# include ecbuild required commands
include(ecbuild_system NO_POLICY_SCOPE)
ecbuild_requires_macro_version(2.0)
ecbuild_declare_project()
ecbuild_enable_fortran(REQUIRED MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module)

include(CompilerFlags)

# build executables
add_subdirectory(src)

# documentation
ecbuild_add_option(FEATURE DOCS
                   DESCRIPTION "Documentation"
                   DEFAULT OFF
                   REQUIRED_PACKAGES "Latex")

if(HAVE_DOCS)
  add_subdirectory(doc)
endif()

# finalize
ecbuild_install_project(NAME dwarf3)

# print summary
ecbuild_print_summary()
