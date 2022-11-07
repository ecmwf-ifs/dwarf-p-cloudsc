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

# set( OpenMP_C_FLAGS "-fopenmp=libomp" CACHE STRING "" )
# set( OpenMP_C_LIB_NAMES "omp;pthread" CACHE STRING "" )
# set( OpenMP_Fortran_FLAGS   "-fopenmp -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx90a" CACHE STRING "" )
# set( OpenMP_Fortran_LIB_NAMES "omp;omptarget;pthread" CACHE STRING "" )
# set( OpenMP_omptarget_LIBRARY "/pfs/lustrep1/projappl/project_465000116/bareuter/rocm/afar/432/llvm/lib/libomptarget.so" CACHE FILEPATH "" )

####################################################################
# OpenAcc FLAGS
####################################################################

set( ENABLE_ACC OFF CACHE STRING "" )

####################################################################
# COMMON FLAGS
####################################################################

set(ECBUILD_Fortran_FLAGS "-v -fpic")
#set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mframe")
#set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mbyteswapio")
#set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mstack_arrays")
#set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mrecursive")
#set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Ktrap=fp")
#set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Kieee")
#set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Mdaz")
