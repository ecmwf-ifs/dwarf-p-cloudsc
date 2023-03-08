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

set( ECBUILD_FIND_MPI OFF CACHE STRING "" )
set( ENABLE_USE_STMT_FUNC ON CACHE STRING "" )

####################################################################
# OpenMP FLAGS
####################################################################

set( ENABLE_OMP ON CACHE STRING "" )
set( OpenMP_C_FLAGS   "-homp" CACHE STRING "" )
set( OpenMP_Fortran_FLAGS   "-homp" CACHE STRING "" )

####################################################################
# OpenACC FLAGS
####################################################################

set( ENABLE_ACC ON CACHE STRING "" )
set( OpenACC_C_FLAGS "-hacc" )
set( OpenACC_CXX_FLAGS "-hacc" )
set( OpenACC_Fortran_FLAGS "-hacc -h acc_model=deep_copy" )

####################################################################
# Compiler FLAGS
####################################################################

# General Flags (add to default)
set(ECBUILD_Fortran_FLAGS "-hcontiguous")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -hbyteswapio")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Wl, --as-needed")

set(ECBUILD_Fortran_FLAGS_BIT "-v")#-O3")# -hfp1 -hscalar3 -hvector3 -haggress -DNDEBUG")

set(HDF5_DIR "/cfs/klemming/home/b/breuter/Public/hdf5/1.12.2/cray/14.0.1" CACHE STRING "")
set(HDF5_ROOT "/cfs/klemming/home/b/breuter/Public/hdf5/1.12.2/cray/14.0.1" CACHE STRING "")
