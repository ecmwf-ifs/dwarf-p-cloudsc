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

set( ECBUILD_FIND_MPI OFF )
set( ENABLE_USE_STMT_FUNC ON CACHE STRING "" )

####################################################################
# OpenMP FLAGS
####################################################################

set( ENABLE_OMP ON CACHE STRING "" )
set( OpenMP_C_FLAGS   "-fopenmp" CACHE STRING "" )
set( OpenMP_CXX_FLAGS "-fopenmp" CACHE STRING "" )
set( OpenMP_Fortran_FLAGS   "-fopenmp -hnoacc -hlist=aimd" CACHE STRING "" )

set( OpenMP_C_LIB_NAMES       "craymp" )
set( OpenMP_CXX_LIB_NAMES     "craymp" )
set( OpenMP_Fortran_LIB_NAMES "craymp" )
set( OpenMP_craymp_LIBRARY    "/opt/cray/pe/cce/15.0.1/cce/x86_64/lib/libcraymp.so" )

####################################################################
# OpenACC FLAGS
####################################################################

set( ENABLE_ACC ON CACHE STRING "" )
set( OpenACC_C_FLAGS "-hacc" )
set( OpenACC_CXX_FLAGS "-hacc" )
set( OpenACC_Fortran_FLAGS "-hacc -h acc_model=deep_copy" )

####################################################################
# OpenACC FLAGS
####################################################################

set(CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} -03 -ffast-math")
if(NOT DEFINED CMAKE_CUDA_ARCHITECTURES)
  set(CMAKE_HIP_ARCHITECTURES gfx90a)
endif()

####################################################################
# Compiler FLAGS
####################################################################

# General Flags (add to default)
set(ECBUILD_Fortran_FLAGS "-hcontiguous")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -hbyteswapio")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Wl, --as-needed")

set(ECBUILD_Fortran_FLAGS_BIT "-O3 -G2 -haggress -DNDEBUG")
# set(ECBUILD_Fortran_FLAGS_BIT "-O3 -hfp1 -hscalar3 -hvector3 -G2 -haggress -DNDEBUG")
