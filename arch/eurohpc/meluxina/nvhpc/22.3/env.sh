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
module_unload ParaStationMPI
module_unload NVHPC
module_unload gompi
module_unload HDF5
module_unload CMake

# Load modules
module use /apps/USE/easybuild/staging/2021.5/modules/all

module_load NVHPC/22.3
module_load ParaStationMPI/5.4.11-1-GCC-10.3.0-CUDA-11.3.1
module_load CMake/3.20.4
module_load CUDA/11.3.1
module_load Boost/1.76.0-GCC-10.3.0
module_load Python/3.9.5-GCCcore-10.3.0

export CC=nvc
export CXX=nvc++
export F77=nvfortran
export FC=nvfortran
export F90=nvfortran

export HDF5_ROOT=/mnt/tier2/project/p200061/nvhpc-install

# Increase stack size to maximum
ulimit -S -s unlimited

set -x

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
