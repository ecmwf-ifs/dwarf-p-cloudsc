This is a driver allowing to execute IFS physics from within a Python script
Steps to run and perform basic test  on ATOS:
Uncomment #add_subdirectory(cloudsc_pyiface) in ./src/CMakeLists.txt (at the moment needed to mitigate CI incompatibility) 
module load python3
python3 -m venv cleanvenv
source ./cleanvenv/bin/activate
./cloudsc-bundle create && ./cloudsc-bundle build --build-type=release --arch=./arch/ecmwf/hpc2020/intel/2021.4.0/ && cd build; . env.sh; ctest --verbose

TODOs:
- should Python be touched in the topmost cmakelists.txt ?
- remove locals
- activate expand routines (currently commented out for no input parameters are available) to enable block memory structure
- investigate why cloudsc loader is needed on Atos but not on Mac
- check performance on Atos

