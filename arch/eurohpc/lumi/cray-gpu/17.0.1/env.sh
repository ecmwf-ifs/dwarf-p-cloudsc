) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Source me to get the correct configure/build/run environment

# Warning: with rocm/6.2.2 there are following liner errors for the executables. With rocm/6.0.3 they are gone
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<long, 1>::assign<true, (int*)0>(long const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<float, 2>::assign<true, (int*)0>(float const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<float, 1>::assign<true, (int*)0>(float const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<double, 1>::assign<true, (int*)0>(double const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<int, 2>::assign<true, (int*)0>(atlas::array::ArrayView<int, 2> const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<float, 3>::assign<true, (int*)0>(float const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<double, 3>::assign<true, (int*)0>(double const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<int, 3>::assign<true, (int*)0>(int const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<int, 1>::assign<true, (int*)0>(atlas::array::ArrayView<int, 1> const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<long, 1>::assign<true, (int*)0>(atlas::array::ArrayView<long, 1> const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<double, 2>::assign<true, (int*)0>(double const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<double, 2>::assign<true, (int*)0>(atlas::array::ArrayView<double, 2> const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<int, 1>::assign<true, (int*)0>(int const&)'


# Store tracing and disable (module is *way* too verbose)
{ tracing_=${-//[^x]/}; set +x; } 2>/dev/null

module_load() {
  echo "+ module load $1"
  if [ "${2:-""}" == "ECBUNDLE_CONFIGURE_ONLY" ]; then
    if [ -n "${ECBUNDLE_CONFIGURE:-""}" ]; then
      module load $1
    else
      echo " WARNING: Module $1 not loaded (only during configuration)"
    fi
  else
    module load $1
  fi
}
module_unload() {
  echo "+ module unload $1"
  module unload $1
}

# Unload to be certain
module reset

# Load modules
module_load LUMI/24.03
module_load partition/G
module_load PrgEnv-cray/8.5.0
module_load cray-mpich/8.1.29
module_load craype-network-ofi
module_load rocm/6.0.3
module_load buildtools/24.03
module_load cray-python/3.10.10
module_load craype-x86-trento
module_load craype-accel-amd-gfx90a
module_load Eigen/3.4.0

### Handling of "magic" cray modules
# 1) Load the cray modules
module_load cray-hdf5/1.12.2.11
# 2) Store variables to locate the packages
_HDF5_ROOT=${CRAY_HDF5_PREFIX}
# 3) Unload the cray modules in reverse order, removing all the magic
module_unload cray-hdf5
# 4) Define variables that CMake introspects
export HDF5_ROOT=${_HDF5_ROOT}

# Export environment variable3s
export MPI_HOME=${MPICH_DIR}
export CC=cc
export CXX=CC
export FC=ftn
export HIPCXX=$(hipconfig --hipclangpath)/clang++

module list

set -x

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
