#!/bin/bash

CC=$(which cc)
CXX=$(which CC)
ENV=nvidia
FORTRAN_MODES=( gpu-scc gpu-scc-hoist )
GT4PY_BACKENDS=( )
HOST=dom
NPROMA=128
NXS=( 512 1024 2048 4096 8192 16384 32768 65536 131072 )

for FORTRAN_MODE in "${FORTRAN_MODES[@]}"
do
  echo "FORTRAN: $FORTRAN_MODE: start"
  for NX in "${NXS[@]}"
  do
    echo -n "  nx=$NX: "
    python run_fortran.py \
      --build-dir=../../../../develop/build/"$ENV" \
      --nproma="$NPROMA" \
      --num-runs=20 \
      --num-threads=1 \
      --output-file=../data/performance_"$ENV".csv \
      --host-alias="$HOST" \
      --mode="$FORTRAN_MODE" \
      --nx="$NX" || true
  done
  echo "FORTRAN: $FORTRAN_MODE: end"
done

for GT4PY_BACKEND in "${GT4PY_BACKENDS[@]}"
do
  echo "Python: $GT4PY_BACKEND: start"
  for NX in "${NXS[@]}"
  do
    echo -n "  nx=$NX: "
    CXX=$CXX CC=$CC python run.py \
      --num-runs=20 \
      --disable-checks \
      --disable-validation \
      --output-file=../data/performance_"$ENV".csv \
      --host-alias="$HOST" \
      --backend="$GT4PY_BACKEND" \
      --nx="$NX" || true
  done
  echo "Python: $GT4PY_BACKEND: end"
done
