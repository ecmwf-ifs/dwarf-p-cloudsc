
[ -d $HOME/arm-sve-tools ] || cp -a ~brx-pridley/arm-sve-tools $HOME
source $HOME/arm-sve-tools/isambard-gcc.bashrc

# Requires manual HDF5 installation, e.g. via the script provided
# in https://git.ecmwf.int/users/nabr/repos/build-scripts 
source $HOME/deps/gcc/11.0.0/deps_env.sh
export HDF5_DIR=$HOME/deps/gcc/11.0.0

export C_INCLUDES=-I$HDF5_DIR/include
export FC_INCLUDES=-I$HDF5_DIR/include

export CC=gcc
export CXX=g++
export FC=gfortran

module load cmake
module load openmpi/4.1.0/gcc-11.0

module use /software/aarch64/tools/arm-compiler/21.0/modulefiles
module load arm21
module load armpl-AArch64-SVE/21.0.0

module list

ulimit -S -s unlimited

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
