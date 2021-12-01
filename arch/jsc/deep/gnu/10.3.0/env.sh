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
module_unload Java
module_unload Python
module_unload HDF5
module_unload Boost
module_unload ParaStationMPI
module_unload GCC
module_unload CMake

# Load modules
module_load CMake/3.18.0
module_load GCC/10.3.0
module_load ParaStationMPI/5.4.9-1
module_load Boost/1.74.0
module_load HDF5/1.10.6
module_load Python/3.8.5
module_load Java/15.0.1

# Increase stack size to maximum
ulimit -S -s unlimited

# Restore tracing to stored setting
if [[ -n "$tracing_" ]]; then set -x; else set +x; fi

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
