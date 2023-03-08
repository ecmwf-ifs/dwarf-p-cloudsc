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

# set( OpenMP_Fortran_FLAGS   "-fopenmp -target x86_64-pc-linux-gnu -fopenmp-targets=amdgcn-amd-amdhsa -Xopenmp-target=amdgcn-amd-amdhsa -march=gfx90a" CACHE STRING "" )
# set( OpenMP_Fortran_FLAGS   "-fopenmp --offload-arch=gfx90a " CACHE STRING "" )
set( OpenMP_Fortran_FLAGS   "-fopenmp --offload-arch=gfx90a" CACHE STRING "" )
# set( OpenMP_OFFLOAD_Fortran_FLAGS "--offload-arch=gfx90a" CACHE STRING "" )

####################################################################
# OpenAcc FLAGS
####################################################################

set( ENABLE_ACC OFF CACHE STRING "" )

####################################################################
# COMMON FLAGS
####################################################################

set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fpic -O3 ${OpenMP_Fortran_FLAGS}")

set(HDF5_DIR "/cfs/klemming/home/b/breuter/Public/hdf5/1.12.2/rocm/5.3.3" CACHE STRING "")
set(HDF5_ROOT "/cfs/klemming/home/b/breuter/Public/hdf5/1.12.2/rocm/5.3.3" CACHE STRING "")
