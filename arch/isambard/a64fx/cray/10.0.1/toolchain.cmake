# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

####################################################################
# Compiler FLAGS
####################################################################

# General Flags (add to default)
set( ECBUILD_FIND_MPI ON )
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -hcontiguous")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -hbyteswapio")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Ktrap=fp")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Wl, --as-needed")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -target-accel=host")

set(ECBUILD_Fortran_FLAGS_BIT "-O3 -hfp1 -hscalar3 -hvector3 -G2 -haggress -DNDEBUG")

####################################################################
# LINK FLAGS
####################################################################

# Compatibility with HDF5 1.12
set(H5_USE_110_API ON)
