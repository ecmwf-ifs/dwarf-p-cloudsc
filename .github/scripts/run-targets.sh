#!/bin/bash
set -euo pipefail
set -x

exit_code=0
cd build

#
# Run each of the binaries using default arguments and validate exit codes
#

for target in $(ls bin)
do
  if [[ "$mpi_flag" == "--with-mpi" ]]
  then
    # Two ranks with one thread each
    mpirun -np 2 bin/$target 1
  else
    # Default arguments
    bin/$target
  fi
  exit_code=$((exit_code + $?))
done

exit $exit_code
