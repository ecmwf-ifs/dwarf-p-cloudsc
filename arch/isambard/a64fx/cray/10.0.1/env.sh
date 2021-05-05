
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
