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
  if [ "$2" == "ECBUILD_CONFIGURE_ONLY" ]; then
    if [ -n "${ECBUILD_CONFIGURE}" ]; then
      echo "+ module load $1"
      module load $1
    else
      echo " WARNING: Module $1 not loaded (only during configuration)"
    fi
  else
    echo "+ module load $1"
    module load $1
  fi
}
module_unload() {
  echo "+ module unload $1"
  module unload $1
}

# Unload to be certain
module_unload cmake
module_unload python
module_unload python3
module_unload boost
module_unload ecbuild
module_unload ifs-support
module_unload cdt
module_unload boost
module_unload PrgEnv-cray
module_unload PrgEnv-intel
module_unload intel
module_unload gcc

export EC_CRAYPE_INTEGRATION=off

prgenvswitchto intel
# Load modules
module_unload intel
module_load intel/2021.1.1
module_load python3
module_load boost
module_load ninja
module_load cmake
module_load boost

set -x

# This is used to download binary test data
export http_proxy="http://slb-proxy-web.ecmwf.int:3333/"

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
