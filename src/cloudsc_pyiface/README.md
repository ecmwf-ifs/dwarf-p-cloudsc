This is a driver allowing to execute IFS physics from within a Python script
Steps to run and perform basic test  on ATOS:
1) Open interactive session on the computing node:
export OMP_NUM_THREADS=64
OMP_PLACES="{$(seq -s '},{' 0 $(($OMP_NUM_THREADS-1)) )}" srun -q np --ntasks=1 --hint=nomultithread --cpus-per-task=$OMP_NUM_THREADS --pty /bin/bash 
module load python3/new
module load pi
python3 -m venv cleanvenv
source ./cleanvenv/bin/activate
./cloudsc-bundle create && ./cloudsc-bundle build --build-type=release --arch=./arch/ecmwf/hpc2020/intel/2021.4.0/ && cd build; ctest --verbose
To test performance, execute:
cd bin/pythonexec
./cloudsc_pydriver.py --numomp=$OMP_NUM_THREADS --ngptot=163840 --nproma=32 
Currently, the performance on ATOS is about ~64400 Mflops/s, which is inferior to the reference result for:
dwarf-cloudsc-fortran-pyiref (~96700) and the dwarf-cloudsc-fortran(~103000)

Similar results can be achieved using GNU compilers on ATOS using module pg and --arch=./arch/ecmwf/hpc2020/gnu/11.2.0/

Troubleshooting:
If the Python script can't find the newly created Python module containing Fortran code, it is very likely that different Python versions were involved during build and execution. Check if correct python module is loaded and if the venv is clean.

TODOs:
- the performance is inferior as compared to the reference ?
- should Python be touched in the topmost cmakelists.txt ?
- why ecbuild_find_python can't find Python libs on ATOS ? Non-critical but surprising.
- find out why expand_r3 does not pass nlev parameter correctly and so expand_r3bis needs to be used
-£ enable compile/run on ATOS with nvhpc. Currently, invalid flags are being passed at the f2py c compilation step, i.e.:
nvc-Error-Unknown switch: -Wno-unused-result
nvc-Error-Unknown switch: -fwrapv
nvc-Error-Unknown switch: -Wno-unused-result
nvc-Error-Unknown switch: -fwrapv
error: Command "nvc -Wno-unused-result -Wsign-compare -DNDEBUG -g -fwrapv -O3 -Wall -O3 -O3 --std=c99 -fPIC -DNPY_DISABLE_OPTIMIZATION=1 -I/home/napz/dwarf-p-cloudsc-test/build/cloudsc-dwarf/src/cloudsc_pyiface/../../module -I/home/napz/dwarf-p-cloudsc-test/build/cloudsc-dwarf/src/cloudsc_pyiface/../common/module/ -I/dev/shm/_tmpdir_.napz.39427566/tmp75szt2xp/src.linux-x86_64-3.10 -I/home/napz/dwarf-p-cloudsc-test/cleanvenv/lib/python3.10/site-packages/numpy/core/include -I/home/napz/dwarf-p-cloudsc-test/cleanvenv/include -I/usr/local/apps/python3/3.10.10-01/include/python3.10 -c /dev/shm/_tmpdir_.napz.39427566/tmp75szt2xp/src.linux-x86_64-3.10/_cloudscmodule.c -o /dev/shm/_tmpdir_.napz.39427566/tmp75szt2xp/dev/shm/_tmpdir_.napz.39427566/tmp75szt2xp/src.linux-x86_64-3.10/_cloudscmodule.o" failed with exit status 1


