#!/bin/bash

python run_fortran.py --num-runs=20 --num-threads=24 --output-file=performance.csv --host-alias=daint --mode=fortran --nx=512    || true
python run_fortran.py --num-runs=20 --num-threads=24 --output-file=performance.csv --host-alias=daint --mode=fortran --nx=1024   || true
python run_fortran.py --num-runs=20 --num-threads=24 --output-file=performance.csv --host-alias=daint --mode=fortran --nx=2048   || true
python run_fortran.py --num-runs=20 --num-threads=24 --output-file=performance.csv --host-alias=daint --mode=fortran --nx=4096   || true
python run_fortran.py --num-runs=20 --num-threads=24 --output-file=performance.csv --host-alias=daint --mode=fortran --nx=8192   || true
python run_fortran.py --num-runs=20 --num-threads=24 --output-file=performance.csv --host-alias=daint --mode=fortran --nx=16384  || true
python run_fortran.py --num-runs=20 --num-threads=24 --output-file=performance.csv --host-alias=daint --mode=fortran --nx=32768  || true
python run_fortran.py --num-runs=20 --num-threads=24 --output-file=performance.csv --host-alias=daint --mode=fortran --nx=65536  || true
python run_fortran.py --num-runs=20 --num-threads=24 --output-file=performance.csv --host-alias=daint --mode=fortran --nx=131072 || true

python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_kfirst --nx=512    || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_kfirst --nx=1024   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_kfirst --nx=2048   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_kfirst --nx=4096   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_kfirst --nx=8192   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_kfirst --nx=16384  || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_kfirst --nx=32768  || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_kfirst --nx=65536  || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_kfirst --nx=131072 || true

python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_ifirst --nx=512    || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_ifirst --nx=1024   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_ifirst --nx=2048   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_ifirst --nx=4096   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_ifirst --nx=8192   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_ifirst --nx=16384  || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_ifirst --nx=32768  || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_ifirst --nx=65536  || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:cpu_ifirst --nx=131072 || true

python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:gpu --nx=512    || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:gpu --nx=1024   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:gpu --nx=2048   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:gpu --nx=4096   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:gpu --nx=8192   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:gpu --nx=16384  || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:gpu --nx=32768  || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:gpu --nx=65536  || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=gt:gpu --nx=131072 || true

python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=cuda --nx=512    || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=cuda --nx=1024   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=cuda --nx=2048   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=cuda --nx=4096   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=cuda --nx=8192   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=cuda --nx=16384  || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=cuda --nx=32768  || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=cuda --nx=65536  || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=cuda --nx=131072 || true

python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=dace:gpu --nx=512    || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=dace:gpu --nx=1024   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=dace:gpu --nx=2048   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=dace:gpu --nx=4096   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=dace:gpu --nx=8192   || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=dace:gpu --nx=16384  || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=dace:gpu --nx=32768  || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=dace:gpu --nx=65536  || true
python run.py --num-runs=20 --disable-checks --disable-validation --output-file=performance.csv --host-alias=daint --backend=dace:gpu --nx=131072 || true
