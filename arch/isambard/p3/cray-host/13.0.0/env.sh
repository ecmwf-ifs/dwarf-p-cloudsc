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

module reset

module use /lustre/projects/bristol/modules/modulefiles

module load PrgEnv-cray/8.2.0
module load craype-accel-host
module load craype-x86-milan
module load cray-pals
module load cray-hdf5/1.12.0.7
module load cmake/3.23.2

module list

ulimit -S -s unlimited

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
