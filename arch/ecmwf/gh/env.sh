#! /usr/bin/env bash
# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Source me to get the correct configure/build/run environment

# Store tracing and disable (module is *way* too verbose)
{ tracing_=${-//[^x]/}; set +x; } 2>/dev/null

module_load() {
  echo "+ module load $1"
  module load $1
}
module_unload() {
  echo "+ module unload $1"
  module unload $1
}

# Unload all modules to be certain
module purge

# Load modules
module use /opt/eb_packages/modules/all

#module_load CMake/3.27.6-NVHPC-24.1-CUDA-12.3.0
module_load Python/3.11.3-GCCcore-12.3.0
module_load OpenMPI/4.1.6-GCC-13.2.0
module_load NVHPC/24.1-CUDA-12.3.0
module_load HDF5/1.14.3-NVHPC-24.1-CUDA-12.3.0
module_load Boost/1.82.0-GCC-12.3.0

export PATH=$PERM/software/gh/cmake/3.30.3/bin:$PATH

# module_load CMake/3.27.6-GCCcore-13.2.0
# module_load NVHPC/24.1-CUDA-12.3.0
# module_load HDF5/1.14.3-NVHPC-24.1-CUDA-12.3.0
# module_load GCC/12.3.0 
export CC=nvc
export CXX=nvc++
#export CUDA_HOST_COMPILER=gcc-11
export CUDAHOSTCXX=nvc++
export FC=nvfortran
export NVHPC_CUDA_HOME=/opt/eb_packages/software/NVHPC/24.1-CUDA-12.3.0/Linux_aarch64/24.1/cuda/12.3

# Increase stack size to maximum
ulimit -S -s unlimited

set -x

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null
