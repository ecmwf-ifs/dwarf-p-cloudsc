#!/usr/bin/env -S bash -lxeo pipefail
cleanup() {
  if [[ -n "${{ PATCH }}" ]]; then
    git reset --hard
    git stash pop
  fi
}

set -u

cd ${{ CLOUDSC_HOME }}

# This syntax is correct because the substitution is being done by JUBE
# and not by bash
if [[ -n "${{ PATCH }}" ]]; then
  git stash
  git apply ${{ PATCH }}
fi

# Trap the exit signal so we always return the repo to the initial state
trap cleanup EXIT

./cloudsc-bundle build --retry-verbose \
    --build-dir=${{ JUBE_WP_ABSPATH }} \
    --arch=${{ ARCH }} \
    ${{ PRECISION_FLAG }} \
    ${{ IO_LIBRARY_FLAG }} \
    ${{ MPI_FLAG }} \
    ${{ OTHER_FLAGS }}
