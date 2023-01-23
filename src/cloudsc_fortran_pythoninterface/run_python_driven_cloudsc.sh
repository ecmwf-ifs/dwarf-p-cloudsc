module load python3
python3 -m venv ~/cstest
source ~/cstest/bin/activate
pip install --upgrade pip
pip install f90wrap h5py 
export DYLD_LIBRARY_PATH=$DYLD_LIBRARY_PATH:.
export LD_LIBRARY_PATH=$LD_LIBRARY_PATH:.
ipython ./dwarfdriver.py
