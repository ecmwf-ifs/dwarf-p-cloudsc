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
module_unload eccodes
module_unload boost
module_unload intel
module_unload cmake
module_unload gnu

# Load modules
module_load intel/18.0.1
module_load boost/1.61.0
module_load cmake/3.15.0

set -x

# Increase stack size to maximum
ulimit -S -s unlimited

# Restore tracing to stored setting
if [[ -n "$tracing_" ]]; then set -x; else set +x; fi

export ECBUILD_TOOLCHAIN=$PWD/../toolchains/ecmwf-leap42-intel.cmake
