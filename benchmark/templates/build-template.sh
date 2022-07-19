#!/usr/bin/env -S bash -xeo pipefail
set -u

cd ${{ CLOUDSC_HOME }}

./cloudsc-bundle build --retry-verbose \
    --build-dir=${{ JUBE_WP_ABSPATH }} \
    --arch=${{ ARCH }} \
    ${{ PRECISION_FLAG }} \
    ${{ IO_LIBRARY_FLAG }} \
    ${{ MPI_FLAG }} \
    ${{ OTHER_FLAGS }}
