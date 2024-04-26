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
module load rhel8/slurm
module load rhel8/global
module load dot
module load dawn-env/2024-04-15
module load intel-oneapi-compilers/2024.1.0
module load intel-oneapi-mpi/2021.12.0
module_load boost/1.84.0
module_load cmake/3.27.9
module_load hdf5/1.14.3

set -x

# below should be the default these days
export EnableImplicitScaling=0
# card 0, tile 0
export ZE_AFFINITY_MASK=0.0
export ONEAPI_DEVICE_SELECTOR=level_zero:gpu
# 256 registers per thread, fewer threads
export SYCL_PROGRAM_COMPILE_OPTIONS="-ze-opt-large-register-file"
# this option affects the overhead of SYCL offload calls, in this case 0 seems to help 
export SYCL_PI_LEVEL_ZERO_USE_IMMEDIATE_COMMANDLISTS=0

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
