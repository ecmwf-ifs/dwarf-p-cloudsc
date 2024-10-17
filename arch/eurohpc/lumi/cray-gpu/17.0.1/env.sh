# Source me to get the correct configure/build/run environment

# Warning: with rocm/6.2.2 there are following liner errors for the executables. With rocm/6.0.3 they are gone
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<long, 1>::assign<true, (int*)0>(long const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<float, 2>::assign<true, (int*)0>(float const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<float, 1>::assign<true, (int*)0>(float const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<double, 1>::assign<true, (int*)0>(double const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<int, 2>::assign<true, (int*)0>(atlas::array::ArrayView<int, 2> const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<float, 3>::assign<true, (int*)0>(float const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<double, 3>::assign<true, (int*)0>(double const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<int, 3>::assign<true, (int*)0>(int const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<int, 1>::assign<true, (int*)0>(atlas::array::ArrayView<int, 1> const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<long, 1>::assign<true, (int*)0>(atlas::array::ArrayView<long, 1> const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<double, 2>::assign<true, (int*)0>(double const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<double, 2>::assign<true, (int*)0>(atlas::array::ArrayView<double, 2> const&)'
#/opt/cray/pe/cce/17.0.1/binutils/x86_64/x86_64-pc-linux-gnu/bin/ld.bfd: ../../../lib/libatlas.so.0.39: undefined reference to `void atlas::array::ArrayView<int, 1>::assign<true, (int*)0>(int const&)'

# Store tracing and disable (module is *way* too verbose)
{ tracing_=${-//[^x]/}; set +x; } 2>/dev/null

module_load() {
  echo "+ module load $1"
  if [ "${2:-""}" == "ECBUNDLE_CONFIGURE_ONLY" ]; then
    if [ -n "${ECBUNDLE_CONFIGURE:-""}" ]; then
      module load $1
    else
      echo " WARNING: Module $1 not loaded (only during configuration)"
    fi
  else
    module load $1
  fi
}
module_unload() {
  echo "+ module unload $1"
  module unload $1
}

# Unload to be certain
module --force purge

# Load modules
module_load CrayEnv
module_load PrgEnv-cray/8.5.0
module_load cce/17.0.1
module_load craype-x86-trento
module_load craype-accel-amd-gfx90a
module_load rocm/6.0.3
module_load cray-mpich/8.1.29
module_load craype-network-ofi
module_load cray-python
module_load Eigen/3.4.0
#module_load cray-dsmml
module_load buildtools/24.03
#export AEC_DIR=/scratch/project_465000454/ifs-deps/libaec-1.1.2

### Handling of "magic" cray modules
# 1) Load the cray modules
module_load cray-libsci/24.03
module_load cray-fftw/3.3.10.7
module_load cray-hdf5
#module_load cray-netcdf
# 2) Store variables to locate the packages
export CRAY_LIBSCI=${CRAY_LIBSCI_PREFIX_DIR}/lib/libsci_cray.so
_FFTW_ROOT=${FFTW_ROOT}
_HDF5_ROOT=${CRAY_HDF5_PREFIX}
#_NETCDF_ROOT=${CRAY_NETCDF_PREFIX}
_MPI_HOME=${MPICH_DIR}
# 3) Unload the cray modules in reverse order, removing all the magic
#module_unload cray-netcdf
module_unload cray-hdf5
module_unload cray-fftw
module_unload cray-libsci
#module_unload cray-mpich
# 4) Define variables that CMake introspects
export FFTW_ROOT=${_FFTW_ROOT}
export HDF5_ROOT=${_HDF5_ROOT}
#export NETCDF_ROOT=${_NETCDF_ROOT}
#export NETCDF_DIR=${_NETCDF_ROOT}

# Export environment variable3s
export MPI_HOME=${_MPI_HOME}
#export CMAKE_TOOLCHAIN_FILE=$PWD/toolchain.cmake
export CC=cc
export CXX=CC
export FC=ftn
export HIPCXX=$(hipconfig --hipclangpath)/clang++
export CRAY_ADD_RPATH=yes
export LIBSCI_ARCH_OVERRIDE=broadwell
  # This is required to work around SIGSEGV in ectrans' SGEMM calls, which
  # occur when "rome" or "milan" are used (backtrace points to openblas_sgemm__naples)

export LDFLAGS="-fuse-ld=bfd"

#. env-rocm.sh
 
export LD_LIBRARY_PATH=$CRAY_LD_LIBRARY_PATH:$LD_LIBRARY_PATH



### Print some exported variables
echo "+ export FFTW_ROOT=${FFTW_ROOT}"
#echo "+ export HDF5_ROOT=${HDF5_ROOT}"
#echo "+ export NETCDF_ROOT=${NETCDF_ROOT}"
echo "+ export MPI_HOME=${MPI_HOME}"
echo "+ export CMAKE_TOOLCHAIN_FILE=${CMAKE_TOOLCHAIN_FILE}"
echo "+ export CC=${CC}"
echo "+ export CXX=${CXX}"
echo "+ export FC=${FC}"
echo "+ export HIPCXX=${HIPCXX}"
echo "+ export CRAY_ADD_RPATH=${CRAY_ADD_RPATH}"
echo "+ export LIBSCI_ARCH_OVERRIDE=${LIBSCI_ARCH_OVERRIDE}"

module list 2>&1

set -x
ulimit -S -s unlimited

# Restore tracing to stored setting
{ if [[ -n "$tracing_" ]]; then set -x; else set +x; fi } 2>/dev/null

