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

set( OpenMP_C_FLAGS             "-fopenmp" )
set( OpenMP_CXX_FLAGS           "-fopenmp" )
set( OpenMP_Fortran_FLAGS       "-fopenmp" )

####################################################################
# COMMON FLAGS
####################################################################

set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fpic")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -flto")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Ofast")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -mcpu=a64fx")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -march=armv8.2-a+sve")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -I$ENV{ARMPL_DIR}/include ")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -DNDEBUG")

set(ECBUILD_C_FLAGS "${ECBUILD_C_FLAGS} -fpic")
set(ECBUILD_CXX_FLAGS "${ECBUILD_CXX_FLAGS} -fpic")
set(ECBUILD_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS} -fpic")
set(ECBUILD_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS} -flto")
set(ECBUILD_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS} -armpl")

set(CMAKE_EXE_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS}")
set(CMAKE_SHARED_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS}")

# Compatibility with HDF5 1.12
set(H5_USE_110_API ON)
