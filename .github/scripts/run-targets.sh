#!/bin/bash
set -euo pipefail
set -x

# These targets don't have an MPI-parallel driver routine
non_mpi_targets=(dwarf-P-cloudMicrophysics-IFSScheme dwarf-cloudsc-c)

exit_code=0
cd build

#
# Run each of the binaries using default arguments and validate exit codes
#

for target in $(ls bin)
do
  if [[ "$mpi_flag" == "--with-mpi" && ! " ${non_mpi_targets[*]} " =~ " $target " ]]
  then
    # Two ranks with one thread each, default NPROMA
    mpirun -np 2 bin/$target 1 100
  else
    # Two threads, default NPROMA
    bin/$target 2 100
  fi
  exit_code=$((exit_code + $?))
done

exit $exit_code
