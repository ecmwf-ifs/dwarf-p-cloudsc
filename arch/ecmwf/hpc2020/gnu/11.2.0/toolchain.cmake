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

####################################################################
# OpenMP FLAGS
####################################################################

# set( OpenMP_Fortran_FLAGS       "-fopenmp -foffload=nvptx-none" )

####################################################################
# OpenAcc FLAGS
####################################################################

# set( OpenACC_Fortran_FLAGS "-fopenacc -foffload=nvptx-none" )

####################################################################
# Compiler FLAGS
####################################################################

# General Flags (add to default)
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -ffpe-trap=invalid,zero,overflow")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fstack-arrays")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fconvert=big-endian")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fbacktrace")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fno-second-underscore")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -ffast-math")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fno-unsafe-math-optimizations")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -march=znver2")

# This is dangerous! But GNU 10+ complains about argument mismatch for MPI routines
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fallow-argument-mismatch")
