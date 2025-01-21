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
  echo "+ module load $*"
  module load $*
}
module_unload() {
  echo "+ module unload $*"
  module unload $*
}
module_purge() {
  echo "+ module purge"
  module purge
}

# Load modules
module_load Stages/2024
module_load OpenSSL/1.1
module_load CUDA/12
module_load StdEnv/2024
module_load NVHPC/24.3-CUDA-12
module_load git/2.41.0-nodocs
module_load CMake/3.26.3
module_load OpenMPI/4.1.6
module_load netCDF/4.9.2
module_load netCDF-Fortran/4.6.1
module_load FFTW/3.3.10
module_load Python/3.11.3
module_load HDF5/1.14.2

# Record the RPATH in the executable
export LD_RUN_PATH=$LD_LIBRARY_PATH

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

