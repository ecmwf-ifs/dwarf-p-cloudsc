#!/bin/bash
set -euo pipefail
set -x

# Set nvhpc version to default value if unset
: "${nvhpc_version:=21.9}"

# Use Atlas' nvhpc installation script
wget https://raw.githubusercontent.com/ecmwf/atlas/develop/tools/install-nvhpc.sh
chmod +x install-nvhpc.sh

# Install nvhpc
./install-nvhpc.sh --version $nvhpc_version --prefix "${GITHUB_WORKSPACE}/nvhpc-install" --tmpdir "${RUNNER_TEMP}"

exit 0
