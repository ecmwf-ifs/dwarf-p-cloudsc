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

export CC=icc
export CXX=icpc
export FC=ifort

module_unload PrgEnv-cray

module use /lustre/projects/bristol/modules/modulefiles/

module load gcc/10.3.0
# module_load intel/oneapi/2021.1
# module_load IntelOneApi/mpi/2021.4.0
source /lustre/software/x86/tools/oneapi-2021.4.0/setvars.sh || true
module_load cmake/3.23.2
export HDF5_DIR=$HOME/dwarf-p-cloudsc/hdf5/intel/2021.4.0
export HDF5_ROOT=$HOME/dwarf-p-cloudsc/hdf5/intel/2021.4.0

set -x

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
