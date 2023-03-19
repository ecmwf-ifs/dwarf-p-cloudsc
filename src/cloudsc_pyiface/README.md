This is a driver allowing to execute IFS physics from within a Python script

Steps to run and perform basic test  on ATOS:
```
# Open interactive session on the computing node:
export OMP_NUM_THREADS=64
OMP_PLACES="{$(seq -s '},{' 0 $(($OMP_NUM_THREADS-1)) )}" srun -q np --ntasks=1 --hint=nomultithread --cpus-per-task=$OMP_NUM_THREADS --pty /bin/bash 

# Build as usual; the PyIface setup will create a custom venv in the build directory
./cloudsc-bundle create
./cloudsc-bundle build --build-type=release --arch=./arch/ecmwf/hpc2020/intel/2021.4.0/

# Go to build and activate the created venv
cd build
. venv_pyiface/bin/activate

# To test performance, execute:
python ./bin/cloudsc_pyiface.py --numomp=$OMP_NUM_THREADS --ngptot=163840 --nproma=32
```
An additional CLI option ``--cloudsc-path=<path-to-python-wrapper>``
can be used if the build location used to run f90wrap has changed.

Currently, the performance on ATOS is about about 64400 Mflops/s, which is inferior to the reference result of:
dwarf-cloudsc-fortran-pyiref (about 96700) and the dwarf-cloudsc-fortran(about 103000)

Similar results can be achieved using GNU compilers on ATOS using module pg and --arch=./arch/ecmwf/hpc2020/gnu/11.2.0/

Troubleshooting:
If the Python script can't find the newly created Python module containing Fortran code, it is very likely that different Python versions were involved during build and execution. Check if correct python module is loaded and if the venv is clean.

TODOs:
- the performance is inferior as compared to the reference ?
- should Python be touched in the topmost cmakelists.txt ? I think current solution is suboptimal. 
- why ecbuild_find_python can't find Python libs on ATOS ? Non-critical but surprising.
- find out why expand_r3 does not pass nlev parameter correctly and so expand_r3bis needs to be used
- enable compile/run on ATOS with nvhpc. Currently, invalid flags are being passed at the f2py c compilation step, i.e.:
nvc-Error-Unknown switch: -Wno-unused-result
nvc-Error-Unknown switch: -fwrapv
nvc-Error-Unknown switch: -Wno-unused-result
nvc-Error-Unknown switch: -fwrapv


