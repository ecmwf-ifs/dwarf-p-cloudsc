#!/usr/bin/env bash

version=1.28.5

TEMPORARY_FILES="${TMPDIR:-/tmp}"
export ECKIT_INSTALL_DIR=$(pwd)/eckit-install
while [ $# != 0 ]; do
    case "$1" in
    "--prefix")
        export ECKIT_INSTALL_DIR="$2"; shift
        ;;
    "--tmpdir")
        TEMPORARY_FILES="$2"; shift
        ;;
    "--version")
        version="$2"; shift
        ;;
    *)
        echo "Unrecognized argument '$1'"
        exit 1
        ;;
    esac
    shift
done

ECKIT_VERSION=${version}
URL=https://github.com/ecmwf/eckit

if [ ! -d "${TEMPORARY_FILES}/eckit" ]; then
  echo "Cloning ${TEMPORARY_FILES}/eckit from [${URL}]"
  mkdir -p ${TEMPORARY_FILES}
  git clone ${URL}
  cd eckit
  git checkout ${ECKIT_VERSION}
else
   echo "Download already present in ${TEMPORARY_FILES}/eckit"
fi

mkdir build
cd build
ecbuild --prefix=$ECKIT_INSTALL_DIR ..
make install