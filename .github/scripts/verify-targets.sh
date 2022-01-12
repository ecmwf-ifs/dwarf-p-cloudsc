#!/bin/bash
set -euo pipefail
set -x

exit_code=0

#
# Build the list of targets
#

targets=(dwarf-P-cloudMicrophysics-IFSScheme)

# Without HDF5/Serialbox, no other variants are currently available
if [[ "$io_library_flag" == "--with-hdf5" ]]
then
  targets+=(dwarf-cloudsc-fortran)

  if [[ "$gpu_flag" == "--with-gpu" ]]
  then
    targets+=(dwarf-cloudsc-gpu-claw dwarf-cloudsc-gpu-scc dwarf-cloudsc-gpu-scc-hoist)
  fi
fi

#
# Verify each target exists
#
echo "::notice title=Targets::Expected targets:\n${targets[@]}"

for target in "${targets[@]}"
do
  if [[ ! -f build/bin/$target ]]
  then
    exit_code=1
    echo "::error::Missing target: $target"
  fi
done

#
# Check there aren't any other binaries
#

if [[ ${#targets[@]} -lt $(ls build/bin | wc -l) ]]
then
  exit_code=1
  echo "::error::Additional targets found in build/bin"
  echo "::error::$(ls build/bin)"
fi

exit $exit_code
