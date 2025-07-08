#!/usr/bin/env bash

version=0.13.3

TEMPORARY_FILES="${TMPDIR:-/tmp}"
export FCKIT_INSTALL_DIR=$(pwd)/fckit-install
while [ $# != 0 ]; do
    case "$1" in
    "--prefix")
        export FCKIT_INSTALL_DIR="$2"; shift
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

FCKIT_VERSION=${version}
URL=https://github.com/ecmwf/fckit

if [ ! -d "${TEMPORARY_FILES}/fckit" ]; then
  echo "Cloning ${TEMPORARY_FILES}/fckit from [${URL}]"
  mkdir -p ${TEMPORARY_FILES}
  git clone ${URL}
  cd fckit
  git checkout ${FCKIT_VERSION}
else
   echo "Download already present in ${TEMPORARY_FILES}/fckit"
fi

mkdir build
cd build
ecbuild --prefix=$FCKIT_INSTALL_DIR ..
make install