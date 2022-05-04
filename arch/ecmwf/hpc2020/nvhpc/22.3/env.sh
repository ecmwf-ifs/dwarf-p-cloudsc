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
module reset

# Load modules
module use /usr/local/apps/hpc_sdk/22.3/modulefiles/
module_load prgenv/nvidia
module_unload nvidia
module_load nvhpc/22.3
source /usr/local/apps/hpc_sdk/22.3/Linux_x86_64/2022/comm_libs/hpcx/latest/hpcx-mt-init.sh hpcx_load
#module_load hdf5/1.10.6
source /home/nabr/hdf5/deps/nvhpc_22.3.sh
export HDF5_DIR=/home/nabr/hdf5/deps/nvhpc/22.3

module_load cmake/3.20.2
module_load python3/3.8.8-01
module_load java/11.0.6

# Increase stack size to maximum
ulimit -S -s unlimited

set -x

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
