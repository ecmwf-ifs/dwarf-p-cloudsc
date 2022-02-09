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
module_unload boost
module_unload cmake
module_unload intel
module_unload pgi
module_unload gnu

# Load modules
module use /opt/nvidia/hpc_sdk/modulefiles
# module load nvhpc
module load nvhpc-nompi/20.9
module_load boost/1.61.0
module_load cmake/3.19.5

set -x

# Increase stack size to maximum
ulimit -S -s unlimited

# Fix boost header location
export BOOST_INCLUDEDIR="/usr/local/apps/boost/1.61.0/PGI/17.1/include/"

# Include local OpenMPI in the path for discovery in build
export PATH="/local/hdd/nabr/openmpi/nvhpc-nompi/20.9/bin:$PATH"

# Custom HDF5 library build with F03 interfaces
export HDF5_ROOT="/local/hdd/nabr/hdf5/nvhpc/20.9"

# Restore tracing to stored setting
if [[ -n "$tracing_" ]]; then set -x; else set +x; fi

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
export ANT_OPTS="-Dhttp.proxyHost=proxy.ecmwf.int -Dhttp.proxyPort=3333 -Dhttps.proxyHost=proxy.ecmwf.int -Dhttps.proxyPort=3333"
