#!/bin/bash
set -euo pipefail
set -x

exit_code=0

#
# Build the list of targets
#

targets=(dwarf-P-cloudMicrophysics-IFSScheme dwarf-cloudsc-fortran dwarf-cloudsc-c)

if [[ "$build_flags" == *"--with-gpu"* ]]
then
  targets+=(dwarf-cloudsc-gpu-scc dwarf-cloudsc-gpu-scc-hoist dwarf-cloudsc-gpu-scc-k-caching)
  targets+=(dwarf-cloudsc-gpu-omp-scc-hoist)
  if [[ "$build_flags" == *"--with-claw"* ]]
  then
    targets+=(dwarf-cloudsc-gpu-claw)
  fi
  if [[ "$build_flags" == *"--with-cuda"* ]]
  then
    targets+=(dwarf-cloudsc-gpu-scc-cuf dwarf-cloudsc-gpu-scc-cuf-k-caching)
    if [[ "$io_library_flag" == "--with-serialbox" ]]
    then
        targets+=(dwarf-cloudsc-c-cuda dwarf-cloudsc-c-cuda-hoist dwarf-cloudsc-c-cuda-k-caching)
    fi
  fi
fi

if [[ "$build_flags" == *"--with-loki"* ]]
then
  targets+=(dwarf-cloudsc-loki-idem dwarf-cloudsc-loki-sca)
  targets+=(dwarf-cloudsc-loki-scc dwarf-cloudsc-loki-scc-hoist)
  targets+=(dwarf-cloudsc-loki-idem-stack dwarf-cloudsc-loki-scc-stack)
  if [[ "$build_flags" != *"--single-precision"* ]]
  then
    targets+=(dwarf-cloudsc-loki-c)
  fi
  if [[ "$build_flags" == *"--with-claw"* ]]
  then
    targets+=(dwarf-cloudsc-loki-claw-cpu dwarf-cloudsc-loki-claw-gpu)
  fi
  if [[ "$build_flags" == *"--with-cuda"* ]]
  then
    targets+=(dwarf-cloudsc-loki-scc-cuf-hoist dwarf-cloudsc-loki-scc-cuf-parametrise)
  fi
fi

if [[ "$build_flags" == *"--with-atlas"* ]]
then
  targets+=(dwarf-cloudsc-fortran-atlas)
fi

if [[ "$build_flags" == *"--cloudsc-fortran-pyiface=ON"* ]]
then
  targets+=(cloudsc_pyiface.py)
fi

if [[ "$build_flags" == *"--cloudsc-python-f2py=ON"* ]]
then
  targets+=(cloudsc_f2py.py)
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

exit $exit_code
