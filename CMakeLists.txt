###############################################################################
# (C) Copyright 1996-2015 ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation nor
# does it submit to any jurisdiction.
################################################################################

cmake_minimum_required( VERSION 2.8.4 FATAL_ERROR )

project( dwarf3-cloud_scheme )

set(ENABLE_OS_TYPES_TEST       OFF   CACHE BOOL " ")
set(ENABLE_OS_FUNCTIONS_TEST   OFF   CACHE BOOL " ")
set(ENABLE_LARGE_FILE_SUPPORT  OFF   CACHE BOOL " ")
set(ENABLE_OS_ENDINESS_TEST    OFF   CACHE BOOL " ")

include(ecbuild_system NO_POLICY_SCOPE)
ecbuild_requires_macro_version(1.9)
ecbuild_declare_project()
ecbuild_enable_fortran(REQUIRED MODULE_DIRECTORY ${CMAKE_BINARY_DIR}/module)

#ecbuild_add_option(FEATURE ATLAS
#                   DESCRIPTION "dwarf2: Support for atlas"
#                   REQUIRED_PACKAGES "PROJECT atlas")
#if(NOT DWARF3_HAVE_ATLAS)
#    ecbuild_warn( "The dwarf3 executables requiring atlas cannot be built")
#endif()

#ecbuild_add_option(FEATURE GRIB
#                   DESCRIPTION "dwarf2: Support for grib data format"
#                   REQUIRED_PACKAGES "PROJECT grib_api VERSION 1.14.3")
#if(NOT DWARF2_HAVE_GRIB)
#    ecbuild_error("The dwarf2 executables require grib_api")
#endif()



###############################################################################
# build executables
##############################################################################
add_subdirectory(src) 



###############################################################################
# finalize
###############################################################################
ecbuild_install_project(NAME dwarf3)
ecbuild_print_summary()
