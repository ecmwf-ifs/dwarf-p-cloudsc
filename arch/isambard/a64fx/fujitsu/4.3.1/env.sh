# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.


[ -d $HOME/arm-sve-tools ] || cp -a ~brx-pridley/arm-sve-tools $HOME
source $HOME/arm-sve-tools/isambard-fujitsu.bashrc

# Requires manual HDF5 installation, e.g. via the script provided
# in https://git.ecmwf.int/users/nabr/repos/build-scripts
source $HOME/deps/fujitsu/4.3.1/deps_env.sh
export HDF5_DIR=$HOME/deps/fujitsu/4.3.1

export CC=fcc
export CXX=FCC
export FC=frt

module load cmake

module list

# export FLIB_HPCFUNC_INFO=TRUE
export FLIB_HPCFUNC=TRUE
export FLIB_FASTOMP=TRUE
# export XOS_MMM_L_PRINT_ENV=on
export XOS_MMM_L_HPAGE_TYPE=hugetlbfs
export XOS_MMM_L_MAX_ARENA_NUM=1

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
