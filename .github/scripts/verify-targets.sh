#!/bin/bash
set -euo pipefail
set -x

exit_code=0

#
# Build the list of targets
#

targets=(dwarf-P-cloudMicrophysics-IFSScheme dwarf-cloudsc-fortran)

if [[ "$io_library_flag" == "--with-serialbox" ]]
then
  targets+=(dwarf-cloudsc-c)
fi

if [[ "$gpu_flag" == "--with-gpu" ]]
then
  targets+=(dwarf-cloudsc-gpu-scc dwarf-cloudsc-gpu-scc-hoist dwarf-cloudsc-gpu-omp-scc-hoist)
  if [[ "$claw_flag" == "--with-claw" ]]
  then
    targets+=(dwarf-cloudsc-gpu-claw)
  fi
fi

if [[ "$loki_flag" == "--with-loki" ]]
then
  targets+=(dwarf-cloudsc-loki-idem dwarf-cloudsc-loki-sca)
  targets+=(dwarf-cloudsc-loki-scc dwarf-cloudsc-loki-scc-hoist)
  if [[ "$prec_flag" != "--single-precision" ]]
  then
    targets+=(dwarf-cloudsc-loki-c)
  fi
  if [[ "$claw_flag" == "--with-claw" ]]
  then
    targets+=(dwarf-cloudsc-loki-claw-cpu dwarf-cloudsc-loki-claw-gpu)
  fi
fi

#
# Verify each target exists
#
echo "::debug::Expected targets: ${targets[@]}"

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
  echo "::error::Expected targets: ${targets[@]}"
  echo "::error::Found targets: $(ls -1 build/bin | tr '\n' ' ')"
fi

exit $exit_code
