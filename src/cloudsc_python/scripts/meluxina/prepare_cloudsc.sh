#!/bin/bash -l

# load env
module load env/release/2022.1

# load compilers and libraries
module load GCC
module load NVHPC/22.7-CUDA-11.7.0
module load CUDA/11.7.0
module load AOCC
module load intel
module load intel-compilers

# load tools
module load Boost
module load CMake
module load HDF5
module load Python

# set local
export LC_ALL=C.UTF-8
export LANG=C.UTF-8

# set/fix CUDA-related variables
NVCC_PATH=$(which nvcc)
CUDA_PATH=$(echo $NVCC_PATH | sed -e "s/\/bin\/nvcc//g")
export CUDA_HOME=$CUDA_PATH
export NVHPC_CUDA_HOME=$CUDA_PATH
export NVHPC_CUDA_LIB_PATH=/apps/USE/easybuild/release/2022.1/software/NVHPC/22.7-CUDA-11.7.0/Linux_x86_64/22.7/compilers/lib
export LD_LIBRARY_PATH="${LD_LIBRARY_PATH/${NVHPC_CUDA_LIB_PATH}:/}"

# required to run OpenACC version on large domains
export PGI_ACC_CUDA_HEAPSIZE=8GB
