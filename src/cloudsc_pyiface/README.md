This is a driver allowing to execute IFS physics from within a Python script, currently adapted to CLOUDSC.

Steps to run and perform basic test  on ATOS:
# Build as usual; the PyIface setup will create a custom venv in the build directory
```
./cloudsc-bundle create
./cloudsc-bundle build --build-type=release  --cloudsc-fortran-pyiface=ON --arch=./arch/ecmwf/hpc2020/intel/2021.4.0/
```
# Work in an interactive session on the computing node:
```
cd build && . env.sh
export OMP_NUM_THREADS=64
OMP_PLACES=cores srun -q np --ntasks=1 --hint=nomultithread --cpus-per-task=$OMP_NUM_THREADS --pty /bin/bash 
```
To test performance, execute:
```
cd build && . env.sh
export OMP_NUM_THREADS=64
./bin/cloudsc_pyiface.py --numomp=$OMP_NUM_THREADS --ngptot=163840 --nproma=32
```
#or, alternatively, submit the non-interactive test job using:
```
OMP_PLACES=cores srun -q np --ntasks=1 --hint=nomultithread --cpus-per-task=$OMP_NUM_THREADS ./bin/cloudsc_pyiface.py --numomp 64 --ngptot 163840 --nproma 32
```

# Additional options
An additional CLI option ``--cloudsc-path=<path-to-python-wrapper>``
can be used if the build location used to run f90wrap has changed.

In addition, to test the Fortran part of the pyiface code independently of the Python driver, 
`` --cloudsc-fortran-pyiface-binary=ON`` option can be used to build Fortran-only binary, mimicking
regular cloudsc-fortran structure. This in particular allows to test if the slight modifications 
to Fortran code alter the computational performance.

# Current performance
Currently, the performance on a single socket with AMD Rome 7742 is about about 64400 Mflops/s, 
which is inferior to the reference result of:
`dwarf-cloudsc-fortran-pyiref` (about 100500), and
`dwarf-cloudsc-fortran` (about 104000)

Similar results can be achieved using GNU compilers on ATOS using `--arch=./arch/ecmwf/hpc2020/gnu/11.2.0/`

# Known issues

### Performance limitations
The performance of PyIface wrapper is inferior as compared to the
`dwarf-cloudsc-fortran` reference. This is probably due to the fact that in
the process of building Fortran binaries, f2py adds low optimization
flags behind the scenes (flags vary between compilers). To
circumevent the problem, a separate explicit compilation step of
f90wrap output files is probably deserved.

### Nvidia compilation
For the same reason, extra effort is needed to enable compile/run on
ATOS with nvhpc. Currently, invalid flags are being passed at the
f2py c compilation step, i.e.:
```
nvc-Error-Unknown switch: -Wno-unused-result
nvc-Error-Unknown switch: -fwrapv
nvc-Error-Unknown switch: -Wno-unused-result
nvc-Error-Unknown switch: -fwrapv
```
