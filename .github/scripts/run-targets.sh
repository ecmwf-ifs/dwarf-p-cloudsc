#!/bin/bash
set -euo pipefail
set -x

# These targets don't have an MPI-parallel driver routine
non_mpi_targets=(dwarf-P-cloudMicrophysics-IFSScheme dwarf-cloudsc-c)

# These targets currently cause issues and are therefore not tested
skipped_targets=(dwarf-cloudsc-gpu-claw)

exit_code=0
cd build

#
# Run each of the binaries with default NPROMA and validate exit codes
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
    # Two ranks with one thread each, default NPROMA
    mpirun -np 2 bin/$target 1 100
  else
    # Single thread, default NPROMA
    bin/$target 1 100
  fi
  exit_code=$((exit_code + $?))
done

exit $exit_code
