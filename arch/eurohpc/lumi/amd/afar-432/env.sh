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

module use /pfs/lustrep1/projappl/project_465000116/bareuter/rocm/mymodules

# Load modules
module_load LUMI/22.06
module_load PrgEnv-aocc/8.3.3
module_load partition/G
module_load rocm/afar-432
module_load buildtools/22.06
module_load cray-python/3.9.12.1

# Load hand-installed dependencies
export HDF5_DIR=/users/bareuter/hdf5/1.12.1/rocm/afar-432
export HDF5_ROOT=$HDF5_DIR

# Specify compilers
export CC=amdclang CXX=amdclang++ FC=amdflang

set -x

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
