#!/usr/bin/env -S bash -lxeo pipefail
set -u
cd build

source ./env.sh

${{ LAUNCH_CMD }} bin/${{ TARGET }} ${{ NUMOMP }} ${{ NGPTOTG }} ${{ NPROMA }} ${{ LAUNCH_CMD_END }}
