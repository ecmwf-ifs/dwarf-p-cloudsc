export LC_ALL=C 
#python3 -m venv ~/cstest
#source ~/cstest/bin/activate
#pip install f90wrap h5py
mkdir build
cd build && rm -rf *
cmake -G Ninja .. && ninja
ninja
cd bin
cd pythonexec
zsh ./run_python_driven_cloudsc.sh
#export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:.
#export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.
#python ./dwarfdriver.py

