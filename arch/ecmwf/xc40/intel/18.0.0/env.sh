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
module_unload grib_api
module_unload eccodes
module_unload emos
module_unload cmake
module_unload fftw
module_unload fcm
module_unload netcdf4
module_unload cray-netcdf-hdf5parallel
module_unload cray-hdf5
module_unload cray-hdf5-parallel
module_unload python
module_unload python3
module_unload boost
module_unload ecbuild
module_unload ifs-support
module_unload fcm
module_unload cdt
module_unload cmake
module_unload boost
module_unload PrgEnv-cray
module_unload PrgEnv-intel
module_unload intel
module_unload gcc

export EC_CRAYPE_INTEGRATION=off

# Load modules
module load gcc
module_load PrgEnv-intel/5.2.82
module_unload intel
module_load intel/18.0.0.033
#module_load gcc/6.3.0
export FFTW_PATH=$( { module load fftw/3.3.4.5; } 2>/dev/null; echo $FFTW_DIR )
export FFTW_LIB_DIR=$FFTW_PATH # for modulefile.in (emos)
export NETCDF_PATH=$( { module load cray-netcdf-hdf5parallel/4.3.2; } 2>/dev/null; echo $NETCDF_DIR )
export HDF5_PATH=$( { module load cray-hdf5-parallel/1.8.13; } 2>/dev/null; echo $HDF5_DIR )
module_load fcm/2015.02.0
module_load python/2.7.12-01
module_load python3/3.6.8-01
module_load boost/1.61.0
module_load eigen/3.2.0
module_load nag
module_load parmetis
module_load cray-snplauncher
module_load atp
module_load ninja
module_load cmake/3.15.0
module_load boost/1.61.0

set -x

# This is used to download binary test data
export http_proxy="http://slb-proxy-web.ecmwf.int:3333/"

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
