# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

export CC=cc
export CXX=CC
export FC=ftn

module use /lustre/projects/bristol/modules/modulefiles

module load PrgEnv-cray/8.0.0
module unload craype-broadwell
module unload craype-network-infiniband
module load craype-accel-host
module load craype-x86-rome
module load cmake/3.23.2
module load intel/mpi/64/2020

export HDF5_DIR=$HOME/dwarf-p-cloudsc/hdf5/cray/11.0.4
export HDF5_ROOT=$HOME/dwarf-p-cloudsc/hdf5/cray/11.0.4

module list

ulimit -S -s unlimited

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
