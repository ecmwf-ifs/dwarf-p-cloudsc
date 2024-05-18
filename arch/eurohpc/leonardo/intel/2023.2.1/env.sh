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

# Load modules
module_load intel-oneapi-compilers/2023.2.1
module_load intel-oneapi-mpi/2021.10.0
module_load cmake/3.24.3
module_load hdf5/1.14.3--intel-oneapi-mpi--2021.10.0--oneapi--2023.2.0
module_load python/3.10.8--gcc--8.5.0

export CC=icc
export CXX=icpc
export F77=ifort
export FC=ifort
export F90=ifort

# Increase stack size to maximum
ulimit -S -s unlimited

set -x

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

# Variable no longer required, make sure it is not set
unset ECBUILD_TOOLCHAIN
