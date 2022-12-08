#!/bin/bash
set -euo pipefail
set -x

hdf5_version=1.10.8

# Choose hdf5
version_parts=($(echo ${hdf5_version} | tr "." "\n"))
major_version=${version_parts[0]}.${version_parts[1]}
temporary_files="${RUNNER_TEMP}/hdf5"
url=https://support.hdfgroup.org/ftp/HDF5/releases/hdf5-${major_version}/hdf5-${hdf5_version}/src/hdf5-${hdf5_version}.tar.gz

# Download hdf5
mkdir -p "${temporary_files}"
curl --location "$url" | tar zx -C "${temporary_files}"

# Build hdf5
cd "${temporary_files}/hdf5-${hdf5_version}"
prefix="${GITHUB_WORKSPACE}/hdf5-install"
mkdir -p "${prefix}"
./configure --prefix="${prefix}" --enable-shared --enable-fortran --enable-hl
make -j
make install

exit 0

