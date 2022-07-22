# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


[ -d $HOME/arm-sve-tools ] || cp -a ~brx-pridley/arm-sve-tools $HOME
source $HOME/arm-sve-tools/isambard-gcc.bashrc

# Requires manual HDF5 installation, e.g. via the script provided
# in https://git.ecmwf.int/users/nabr/repos/build-scripts
source $HOME/deps/gcc/11.0.0/deps_env.sh
export HDF5_DIR=$HOME/deps/arm/21.0.0

export C_INCLUDES=-I$HDF5_DIR/include
export FC_INCLUDES=-I$HDF5_DIR/include

export CC=gcc
export CXX=g++
export FC=gfortran

module load cmake
module load openmpi/4.1.0/gcc-11.0

# ARM PL math library
module use /software/aarch64/tools/arm-compiler/21.0/modulefiles
module load arm21
module load armpl-AArch64-SVE/21.0.0

module list

ulimit -S -s unlimited

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
