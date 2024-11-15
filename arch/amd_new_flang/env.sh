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
# module reset

# Load modules
module_load LUMI/23.09
module_load partition/G
# module_load PrgEnv-cray/8.4.0
module_load cce/16.0.1
# module_load cray-mpich/8.1.27
# module_load craype-network-ofi
# module_load rocm/5.2.3
# module_load buildtools/23.09
# module_load Boost/1.82.0-cpeCray-23.09
# module_load cray-python/3.10.10
module_load craype-x86-trento
module_load craype-accel-amd-gfx90a

### Handling of "magic" cray modules
# # 1) Load the cray modules
# module_load cray-hdf5/1.12.2.7
# # 2) Store variables to locate the packages
# _HDF5_ROOT=${CRAY_HDF5_PREFIX}
# # 3) Unload the cray modules in reverse order, removing all the magic
# module_unload cray-hdf5
# # 4) Define variables that CMake introspects
# export HDF5_ROOT=${_HDF5_ROOT}

# . /scratch/project_465000527/staneker/rocm-afar-6356-drop-4.1.0/env.sh

# spack load aocc
# spack load hdf5

export HDF5_ROOT=/scratch/project_465000527/staneker/flang-new
export HDF5_DIR=/scratch/project_465000527/staneker/flang-new
export HDF5_INCLUDE_DIRS=/scratch/project_465000527/staneker/flang-new/include
export HDF5_LIBS=/scratch/project_465000527/staneker/flang-new/lib

# export HDF5_ROOT=/pfs/lustrep4/scratch/project_465000527/staneker/spack/opt/spack/linux-sles15-zen2/aocc-5.0.0/hdf5-1.14.5-icewtbjjzt3l46whl5nwwxp36x2xaxao
# export HDF5_DIR=/pfs/lustrep4/scratch/project_465000527/staneker/spack/opt/spack/linux-sles15-zen2/aocc-5.0.0/hdf5-1.14.5-icewtbjjzt3l46whl5nwwxp36x2xaxao/cmake
# export HDF5_INCLUDE_DIRS=/pfs/lustrep4/scratch/project_465000527/staneker/spack/opt/spack/linux-sles15-zen2/aocc-5.0.0/hdf5-1.14.5-icewtbjjzt3l46whl5nwwxp36x2xaxao/include
# export HDF5_LIBS=/pfs/lustrep4/scratch/project_465000527/staneker/spack/opt/spack/linux-sles15-zen2/aocc-5.0.0/hdf5-1.14.5-icewtbjjzt3l46whl5nwwxp36x2xaxao/lib
# export HDF5_C_IS_PARALLEL=0

# export OpenMP_ROOT=/pfs/lustrep4/scratch/project_465000527/staneker/spack/opt/spack/linux-sles15-zen2/gcc-13.2.1/aocc-5.0.0-nppydvagz3joo4xdlrgjte7h6qmy2fna/lib

. /scratch/project_465000527/staneker/rocm-afar-6356-drop-4.1.0/env.sh

# Export environment variable3s
# export MPI_HOME=${MPICH_DIR}
# export cc=clang
# export CC=clang++
# export ftn=flang
# export CC=cc
# export CXX=CC
# export FC=ftn
export CC=amdclang CXX=amdclang++ FC=amdflang
# export HIPCXX=$(hipconfig --hipclangpath)/clang++

# module list

# set -x
# set +x

# Restore tracing to stored setting
# { if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
