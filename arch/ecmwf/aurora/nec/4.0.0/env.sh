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

export FC=nfort
export CC=ncc
export CXX=nc++

set -x

# Increase stack size to maximum
ulimit -S -s unlimited

# Enable floating point error trapping at run time
export VE_FPE_ENABLE=DIV,INV,FOF,FUF,INE


export PATH="/local/hdd/nabr/openmpi/nvhpc-nompi/20.9/bin:$PATH"

# Restore tracing to stored setting
if [[ -n "$tracing_" ]]; then set -x; else set +x; fi

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
