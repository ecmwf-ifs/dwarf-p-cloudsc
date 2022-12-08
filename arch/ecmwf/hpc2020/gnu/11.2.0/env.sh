# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

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
module_unload gcc
module_unload openmpi
module_unload hpcx-openmpi
module_unload boost
module_unload hdf5
module_unload cmake
module_unload python3
module_unload java

# Load modules
module_load prgenv/gnu
module_load gcc/11.2.0
module_load hpcx-openmpi/2.10.0
module_load boost/1.71.0
module_load hdf5/1.10.6
module_load cmake/3.20.2
module_load python3/3.8.8-01
module_load java/11.0.6

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
