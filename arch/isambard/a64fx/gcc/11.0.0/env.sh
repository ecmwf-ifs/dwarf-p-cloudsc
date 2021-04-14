
[ -d $HOME/arm-sve-tools ] || cp -a ~brx-pridley/arm-sve-tools $HOME
source $HOME/arm-sve-tools/isambard-gcc.bashrc

module swap gcc/11-20210321
export CC=gcc
export CXX=g++
export FC=gfortran

module load cmake
module load hdf5/1.12.0/gcc-11
module load openmpi/4.1.0/gcc-11.0

module use /software/aarch64/tools/arm-compiler/21.0/modulefiles
module load arm21
module load armpl-AArch64-SVE/21.0.0

module list

ulimit -S -s unlimited

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
