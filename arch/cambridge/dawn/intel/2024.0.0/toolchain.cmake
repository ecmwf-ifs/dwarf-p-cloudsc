# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

####################################################################
# COMPILER FLAGS
####################################################################

set( OpenMP_Fortran_FLAGS "-fopenmp" CACHE STRING "" )
set( OpenMP_C_FLAGS "-fopenmp" CACHE STRING "" )
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -std=c++17 -O3")

####################################################################
# LINK FLAGS
####################################################################

set( ECBUILD_SHARED_LINKER_FLAGS "-Wl,--eh-frame-hdr " )
set( ECBUILD_MODULE_LINKER_FLAGS "-Wl,--eh-frame-hdr -Wl,-Map,loadmap" )
set( ECBUILD_EXE_LINKER_FLAGS    "-Wl,--eh-frame-hdr -Wl,-Map,loadmap -Wl,--as-needed" )
