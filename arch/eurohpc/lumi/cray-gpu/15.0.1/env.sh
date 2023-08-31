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
module_load PrgEnv-cray/8.3.3
module_load LUMI/23.03
# module_load partition/G
module_load rocm/5.2.3
module_load cce/15.0.1
module_load cray-libsci/22.08.1.1
module_load cray-mpich/8.1.18
module_load craype/2.7.20
module_load craype-accel-amd-gfx90a
module_load buildtools/23.03
module_load cray-hdf5/1.12.1.5
module_load cray-python/3.9.12.1
module_load Boost/1.81.0-cpeCray-23.03
module_load partition/G

module list

set -x

export CC=cc CXX=CC FC=ftn

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
