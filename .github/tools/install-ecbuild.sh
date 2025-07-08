#!/usr/bin/env bash

version=3.9.1

TEMPORARY_FILES="${TMPDIR:-/tmp}"
export ECBUILD_INSTALL_DIR=$(pwd)/ecbuild-install
while [ $# != 0 ]; do
    case "$1" in
    "--prefix")
        export ECBUILD_INSTALL_DIR="$2"; shift
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

ECBUILD_VERSION=${version}
URL=https://github.com/ecmwf/ecbuild

if [ ! -d "${TEMPORARY_FILES}/ecbuild" ]; then
  echo "Cloning ${TEMPORARY_FILES}/ecbuild from [${URL}]"
  mkdir -p ${TEMPORARY_FILES}
  git clone ${URL}
  cd ecbuild
  git checkout ${ECBUILD_VERSION}
else
   echo "Download already present in ${TEMPORARY_FILES}/ecbuild"
fi

cd ${TEMPORARY_FILES}/ecbuild
mkdir -p bootstrap && cd bootstrap && rm -rf *
../bin/ecbuild --prefix=$ECBUILD_INSTALL_DIR ..
make install