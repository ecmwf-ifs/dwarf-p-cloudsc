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
set( ENABLE_USE_STMT_FUNC ON CACHE STRING "" )

####################################################################
# OpenACC FLAGS
####################################################################

set( ENABLE_ACC OFF CACHE STRING "" )
set( OpenACC_C_FLAGS "-hnoacc" )
set( OpenACC_CXX_FLAGS "-hnoacc" )
set( OpenACC_Fortran_FLAGS "-hnoacc" )

####################################################################
# Compiler FLAGS
####################################################################

# General Flags (add to default)
set(ECBUILD_Fortran_FLAGS "-hcontiguous")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -hbyteswapio")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Ktrap=fp")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Wl, --as-needed")

set(ECBUILD_Fortran_FLAGS_BIT "-emf -N 1023 -O3 -hfp1 -hscalar3 -hvector3 -G2 -haggress -DNDEBUG")
