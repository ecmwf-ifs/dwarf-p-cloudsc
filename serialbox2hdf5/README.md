# serialbox2hdf5

Convert [Serialbox](https://gridtools.github.io/serialbox/) data to HDF5.

## Quickstart

### 0. Load required modules (e.g. for ECMWF workstations)

```
module load python3
module load hdf5
module load cmake
```

### 1. Create virtual environment

```
python3 -m venv venv
. venv/bin/activate
```

### 2. Install HDF5

```
pip install --upgrade pip wheel
pip install h5py
```

### 3. Install Serialbox

```
git clone https://github.com/gridtools/serialbox.git
mkdir serialbox/build
cd serialbox/build
cmake ../ -DSERIALBOX_ENABLE_PYTHON=ON -DCMAKE_INSTALL_PREFIX=../../venv
cmake --build .
cmake --build . --target install
cd ../..
```

### 4. Run

```
PYTHONPATH=venv/python:$PYTHONPATH python3 serialbox2hdf5.py
```
