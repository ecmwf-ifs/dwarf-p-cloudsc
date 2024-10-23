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

# Unload to be certain
module reset

# Load modules
# module_load LUMI/22.08
# module_load LUMI/23.09
module_load LUMI/24.03
module_load partition/G
# module_load PrgEnv-aocc/8.3.3
# module_load PrgEnv-amd/8.3.3
module_load PrgEnv-amd/8.5.0
# module_load rocm/5.2.5
# module_load rocm/5.3.3
# module_load rocm/5.4.6
# module_load rocm/5.6.1
module_load rocm/6.2.2
# module_load rocm/6.2.2
module_load craype-accel-amd-gfx90a
# module_load buildtools/22.08
# module_load buildtools/23.09
module_load buildtools/24.03
# module_load cray-hdf5/1.12.1.5
# module_load cray-hdf5/1.12.2.11
# module_load cray-python/3.9.12.1

### Handling of "magic" cray modules
# 1) Load the cray modules
# module_load cray-hdf5/1.12.2.7
module_load cray-hdf5/1.12.2.11
# 2) Store variables to locate the packages
_HDF5_ROOT=${CRAY_HDF5_PREFIX}
# 3) Unload the cray modules in reverse order, removing all the magic
module_unload cray-hdf5
# 4) Define variables that CMake introspects
export HDF5_ROOT=${_HDF5_ROOT}

# Specify compilers
export CC=amdclang CXX=amdclang++ FC=amdflang
# export CC=clang CXX=clang++ FC=flang
# export CC=cc CXX=CC FC=ftn
export HIPCXX=$(hipconfig --hipclangpath)/clang++
# export HIPCXX=$(hipconfig --hipclangpath)/bin/hipcc

set -x

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
