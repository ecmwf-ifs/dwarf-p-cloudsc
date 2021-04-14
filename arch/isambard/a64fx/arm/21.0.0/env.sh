
[ -d $HOME/arm-sve-tools ] || cp -a ~brx-pridley/arm-sve-tools $HOME
source $HOME/arm-sve-tools/isambard-arm.bashrc

# Requires manual HDF5 installation, e.g. via the script provided
# in https://git.ecmwf.int/users/nabr/repos/build-scripts 
source $HOME/deps/arm/21.0.0/deps_env.sh
export HDF5_DIR=$HOME/deps/arm/21.0.0

export C_INCLUDES=-I$HDF5_DIR/include
export FC_INCLUDES=-I$HDF5_DIR/include

export CC=armclang
export CXX=armclang++
export FC=armflang
 
module load openmpi/4.1.0/arm-21.0
module load cmake

export ECBUILD_TOOLCHAIN="./toolchain.cmake"
