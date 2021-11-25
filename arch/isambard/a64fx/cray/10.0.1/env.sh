# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


[ -d $HOME/arm-sve-tools ] || cp -a ~brx-pridley/arm-sve-tools $HOME
source $HOME/arm-sve-tools/isambard-cray.bashrc

export CC=cc
export CXX=CC
export FC=ftn

module use /lustre/projects/bristol/modules-a64fx/modulefiles
module load cmake
module load cray-hdf5/1.12.0.2
module load cray-mvapich2_noslurm_nogpu/2.3.4

module list

ulimit -S -s unlimited

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
