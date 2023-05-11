#!/bin/bash
set -eu
set -x

# These targets don't have an MPI-parallel driver routine
non_mpi_targets=(dwarf-P-cloudMicrophysics-IFSScheme dwarf-cloudsc-c)

# These targets currently cause issues and are therefore not tested
skipped_targets=(dwarf-cloudsc-gpu-claw)

if [[ "$arch" == *"nvhpc"* ]]
then
  # Skip GPU targets if built with nvhpc (don't have GPU in test runner)
  skipped_targets+=(dwarf-cloudsc-gpu-scc dwarf-cloudsc-gpu-scc-hoist dwarf-cloudsc-gpu-scc-k-caching) 
  skipped_targets+=(dwarf-cloudsc-gpu-omp-scc-hoist dwarf-cloudsc-gpu-scc-field)

  # Skip GPU targets from Loki if built with nvhpc (don't have GPU in test runner)
  skipped_targets+=(dwarf-cloudsc-loki-claw-gpu dwarf-cloudsc-loki-scc dwarf-cloudsc-loki-scc-hoist)

  # Skip CUDA targets if built with nvhpc
  skipped_targets+=(dwarf-cloudsc-gpu-scc-cuf dwarf-cloudsc-gpu-scc-cuf-k-caching)
  skipped_targets+=(dwarf-cloudsc-loki-scc-cuf-hoist dwarf-cloudsc-loki-scc-cuf-parametrise)
  skipped_targets+=(dwarf-cloudsc-cuda dwarf-cloudsc-cuda-hoist dwarf-cloudsc-cuda-k-caching)

  # Skip C target if built with nvhpc, segfaults for unknown reasons
  skipped_targets+=(dwarf-cloudsc-c dwarf-cloudsc-loki-c)
fi

exit_code=0
cd build

#
# Run each of the binaries with a safe NPROMA value and validate exit codes
#

for target in $(ls bin)
do
  # Skip some targets
  if [[ " ${skipped_targets[*]} " =~ " $target " ]]
  then
    continue
  fi

  if [[ "$mpi_flag" == "--with-mpi" && ! " ${non_mpi_targets[*]} " =~ " $target " ]]
  then
    # Two ranks with one thread each, safe NPROMA
    # NB: Use oversubscribe to run, even if we end up on a single core agent
    mpirun --oversubscribe -np 2 bin/$target 1 100 64
  elif [[ "$target" == "cloudsc_pyiface.py" ]]
  then
    bin/$target --numomp 1 --ngptot 100 --nproma 64
  elif [[ "$target" == "cloudsc_f2py.py" ]]
  then
    bin/$target --ngptot 100 --nproma 128
  else
    # Single thread, safe NPROMA
    bin/$target 1 100 64
  fi
  exit_code=$((exit_code + $?))
done

exit $exit_code
