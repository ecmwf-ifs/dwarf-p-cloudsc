# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

####################################################################
# OpenMP FLAGS
####################################################################

set( OpenMP_C_FLAGS           "-fopenmp" CACHE STRING "" )
set( OpenMP_CXX_FLAGS         "-fopenmp" CACHE STRING "" )
set( OpenMP_Fortran_FLAGS     "-fopenmp" CACHE STRING "" )

####################################################################
# OpenACC FLAGS
####################################################################

set( OpenACC_C_FLAGS "-hacc" CACHE STRING "" )

####################################################################
# Compiler FLAGS
####################################################################

# General Flags (add to default)
set(ECBUILD_Fortran_FLAGS "${ECBUILD_FORTRAN_FLAGS} -hcontiguous")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -hbyteswapio")

set(ECBUILD_Fortran_FLAGS_BIT "-O3 -hfp1 -hscalar3 -hvector3 -G2 -haggress -DNDEBUG")

set(CMAKE_HIP_FLAGS "-Wno-unused-result")

if(NOT DEFINED CMAKE_HIP_ARCHITECTURES)
  set(CMAKE_HIP_ARCHITECTURES gfx90a)
endif()

# select OpenMP pragma to be used
set( HAVE_OMP_TARGET_LOOP_CONSTRUCT OFF CACHE BOOL "" )
set( HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_PARALLEL OFF CACHE BOOL "" )
set( HAVE_OMP_TARGET_LOOP_CONSTRUCT_BIND_THREAD OFF CACHE BOOL "" )
