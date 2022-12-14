# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

####################################################################
# COMPILER
####################################################################

set( ECBUILD_FIND_MPI ON )
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

####################################################################
# OpenMP FLAGS
####################################################################

set( OpenMP_C_FLAGS   "-fopenmp" CACHE STRING "" )
set( OpenMP_CXX_FLAGS   "-fopenmp" CACHE STRING "" )
set( OpenMP_Fortran_FLAGS   "-fopenmp" CACHE STRING "" )

####################################################################
# COMMON FLAGS
####################################################################

foreach(LANG IN ITEMS C CXX Fortran LINKER)
    set(ECBUILD_${LANG}_FLAGS "${ECBUILD_${LANG}_FLAGS} -fpic")
    set(ECBUILD_${LANG}_FLAGS "${ECBUILD_${LANG}_FLAGS} -flto")
endforeach()

set(ECBUILD_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS} -armpl")
set(CMAKE_EXE_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS}")

set(ECBUILD_Fortran_FLAGS_BIT "-Ofast -mcpu=a64fx -march=armv8.2-a+sve -I$ENV{ARMPL_INCLUDES} -DNDEBUG")

# Compatibility with HDF5 1.12
set(H5_USE_110_API ON)
