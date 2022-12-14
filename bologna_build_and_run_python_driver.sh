export LC_ALL=C 
module purge
module load pg ecbuild ninja hdf5 python3 gcc/11.2.0 
python3 -m venv ~/cstest
source ~/cstest/bin/activate
pip install --upgrade pip
pip install f90wrap h5py
make build
cd build
rm -rf *
cmake -G Ninja .. && ninja
ninja
cd bin
cd pythonexec
sh ./run_python_driven_cloudsc.sh

