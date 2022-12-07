This folder contains a Python implementation of the CLOUDSC microphysics scheme based on
[GT4Py](https://github.com/GridTools/gt4py/tree/master). The code is bundled as an installable
package called `cloudsc4py`, whose source code is placed under `src/`.

We strongly recommend installing the package in an isolated virtual environment, which can be
created by issuing the following command from within this directory:
```shell
$ python -m venv venv
```
The virtual environment will be contained in the folder `venv/` and can be activated with
```shell
$ source venv/bin/activate
```
and deactivated with
```shell
$ (venv) deactivate
```
The package `cloudsc4py` can be installed via the Python package manager [pip](https://pypi.org/project/pip/):
```shell
$ (venv) pip install -e .
```
The resulting installation will work on CPU only. To get access to the GPU-accelerated backends of
GT4Py, [CuPy](https://cupy.dev/) is required. We suggest installing CuPy as a precompiled binary
package (wheel)
```shell
$ (venv) pip install cupy-cudaXXX  # XXX stands for the CUDA version available on the system
```
If the installation of CuPy completed successfully, the command
```shell
$ (venv) python -c "import cupy"
```
should produce no output.
All the aforementioned steps can be executed in a single shot by executing the Bash script `bootstrap_venv.sh`:
```shell
$ FRESH_INSTALL=1 VENV=venv INSTALL_CUPY=1 CUPY_VERSION=cupy-cudaXXX [PIP_UPGRADE=1 INSTALL_PRE_COMMIT=1] ./bootstrap_venv.sh
```

The scheme comes in two forms: one where computations are carried out in a single stencil
(see `src/cloudsc4py/{physics,_stencils}/cloudsc.py`), and one where calculations are split into two
stencils (one computing tendencies on the main vertical levels, the other computing fluxes at the
interface levels; see `src/cloudsc4py/{physics,_stencils}/cloudsc_split.py`).

The easiest way to run the dwarf is through the driver scripts `drivers/run.py` and `drivers/run_split.py`.
Run the two scripts with `--help` option to get the full list of command-line options.

For the sake of convenience, we provide the driver `drivers/run_fortran.py` to invoke one of the
FORTRAN variants of the dwarf from Python, and the Bash script `drivers/run_batch.sh` to run the
FORTRAN and Python implementations under different settings.
