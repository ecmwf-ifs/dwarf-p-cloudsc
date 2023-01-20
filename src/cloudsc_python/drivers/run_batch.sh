# (C) Copyright 2018- ECMWF.
# (C) Copyright 2022- ETH Zurich.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

#!/bin/bash

# === general
# name of the host machine
HOST=meluxina
# list of number of columns
NUM_COLS_L=( 512 1024 2048 4096 8192 16384 32768 65536 131072 262144 )

# === FORTRAN
# list of environments
# options: nvhpc
FORTRAN_ENV_L=( )
# list of variants
# options: fortran, gpu-scc, gpu-scc-hoist, gpu-omp-scc-hoist
FORTRAN_VARIANT_L=( fortran gpu-scc gpu-scc-hoist gpu-omp-scc-hoist )
# list of NPROMA values (array must have the same length of FORTRAN_VARIANT_L)
# recommended values: 32 for CPUs, 128 on GPUs
NPROMA_L=( 32 128 128 128 )
# list of number of threads (array must have the same length of FORTRAN_VARIANT_L)
# recommended values: 24 on Piz Daint's CPUs, 128 on MLux's CPUs, 1 on GPUs
FORTRAN_NUM_THREADS_L=( 128 1 1 1 )

# === python
# list of environments
# options: aocc gcc intel
PYTHON_ENV_L=( aocc gcc intel )
# list of C compilers (array must have the same length of PYTHON_ENV_L)
CC_L=( clang gcc icx )
# list of C++ compilers (array must have the same length of PYTHON_ENV_L)
CXX_L=( clang++ g++ icx )
# list of C++ compiler flags (array must have the same length of PYTHON_ENV_L)
CXXFLAGS_L=( "-fbracket-depth=1024" "" "-fbracket-depth=1024" )
# list of linker flags (array must have the same length of PYTHON_ENV_L)
LFLAGS_L=( "" "" "-lstdc++" )
# list of GT4Py backends
# options: numpy, gt:cpu_ifirst, gt:cpu_kfirst, gt:gpu, cuda, dace:cpu, dace:gpu
GT4PY_BACKEND_L=( gt:cpu_ifirst gt:cpu_kfirst dace:cpu )
# list of number of threads (array must have the same length of GT4PY_BACKEND_L)
# recommended values: 24 on Piz Daint, 128 on MLux
PYTHON_NUM_THREADS_L=( 128 128 128 128 128 )

echo "FORTRAN: start"
LEN_FORTRAN_ENV_L=${#FORTRAN_ENV_L[@]}
LEN_FORTRAN_VARIANT_L=${#FORTRAN_VARIANT_L[@]}

for (( i=0; i<"$LEN_FORTRAN_ENV_L"; i++ )); do
  ENV=${FORTRAN_ENV_L[$i]}
  echo "  Env: $ENV: start"

  for (( j=0; j<"$LEN_FORTRAN_VARIANT_L"; j++ )); do
    VARIANT=${FORTRAN_VARIANT_L[$j]}
    mkdir -p ../data/"$HOST"/"$ENV"
    echo "    Variant: $VARIANT: start"
    for NUM_COLS in "${NUM_COLS_L[@]}"; do
      echo -n "      num_cols=$NUM_COLS: "
      python run_fortran.py \
        --build-dir=../../../../develop/build/"$ENV" \
        --nproma="${NPROMA_L[$j]}" \
        --num-runs=20 \
        --num-threads="${FORTRAN_NUM_THREADS_L[$j]}" \
        --output-csv-file=../data/"$HOST"/"$ENV"/performance.csv \
        --host-alias="$HOST" \
        --variant="$VARIANT" \
        --num-cols="$NUM_COLS" || true
    done
    echo "    Variant: $FORTRAN_MODE: end"
  done
  echo "  Env: $ENV: end"
done
echo "FORTRAN: end"

echo ""

echo "Python: start"
LEN_PYTHON_ENV_L=${#PYTHON_ENV_L[@]}
LEN_GT4PY_BACKEND_L=${#GT4PY_BACKEND_L[@]}

for (( i=0; i<"$LEN_PYTHON_ENV_L"; i++ )); do
  ENV=${PYTHON_ENV_L[$i]}
  echo "  Env: $ENV: start"
  export GT_CACHE_ROOT=$PWD/../gt_cache/"$ENV"
  mkdir -p ../data/"$HOST"/"$ENV"

  for (( j=0; j<"$LEN_GT4PY_BACKEND_L"; j++ )); do
    GT4PY_BACKEND=${GT4PY_BACKEND_L[$j]}
    echo "    Backend: $GT4PY_BACKEND: start"

    for NUM_COLS in "${NUM_COLS_L[@]}"; do
      echo -n "      num_cols=$NUM_COLS: "
      OMP_NUM_THREADS=${PYTHON_NUM_THREADS_L[$j]} \
      CC=${CC_L[$i]} CXX=${CXX_L[$i]} CXXFLAGS=${CXXFLAGS_L[$i]} LFLAGS=${LFLAGS_L[$i]} CUDA_HOST_CXX=${CXX_L[$i]} \
        python run_split.py \
        --num-runs=20 \
        --disable-checks \
        --disable-validation \
        --host-alias="$HOST" \
        --backend="$GT4PY_BACKEND" \
        --num-cols="$NUM_COLS" \
        --output-csv-file=../data/"$HOST"/"$ENV"/performance_split.csv \
        --output-csv-file-stencils=../data/"$HOST"/"$ENV"/performance_split_stencils.csv || true
    done
    echo "    Backend: $GT4PY_BACKEND: end"
  done
  echo "  Env: $ENV: end"
done
echo "Python: end"
