#!/usr/bin/env bash

version=0.43.0

TEMPORARY_FILES="${TMPDIR:-/tmp}"
export ATLAS_INSTALL_DIR=$(pwd)/atlas-install
while [ $# != 0 ]; do
    case "$1" in
    "--prefix")
        export ATLAS_INSTALL_DIR="$2"; shift
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

ATLAS_VERSION=${version}
URL=https://github.com/ecmwf/atlas

if [ ! -d "${TEMPORARY_FILES}/atlas" ]; then
  echo "Cloning ${TEMPORARY_FILES}/atlas from [${URL}]"
  mkdir -p ${TEMPORARY_FILES}
  git clone ${URL}
  cd atlas
  git checkout ${ATLAS_VERSION}
else
   echo "Download already present in ${TEMPORARY_FILES}/atlas"
fi

ecbuild --prefix=$ATLAS_INSTALL_DIR -- ${TEMPORARY_FILES}/atlas

make -j10
make install