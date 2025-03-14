# dwarf-p-cloudsc

[![license](https://img.shields.io/github/license/ecmwf-ifs/dwarf-p-cloudsc)](https://www.apache.org/licenses/LICENSE-2.0.html)
[![build](https://github.com/ecmwf-ifs/dwarf-p-cloudsc/actions/workflows/build.yml/badge.svg)](https://github.com/ecmwf-ifs/dwarf-p-cloudsc/actions/workflows/build.yml)

`dwarf-p-cloudsc` is intended to test the CLOUDSC cloud microphysics scheme of the IFS.

*This package is made available to support research collaborations and is not
officially supported by ECMWF*

## Contact

Michael Lange (michael.lange@ecmwf.int),
Willem Deconinck (willem.deconinck@ecmwf.int),
Balthasar Reuter (balthasar.reuter@ecmwf.int)

## Licence

`dwarf-p-cloudsc` is distributed under the Apache Licence Version 2.0. See
[LICENSE](LICENSE) file for details.

## Contributing

Contributions to `dwarf-p-cloudsc` are welcome. 
In order to do so, please create a pull request with your contribution and sign the [contributors license agreement (CLA)](https://bol-claassistant.ecmwf.int/ecmwf-ifs/dwarf-p-cloudsc).

## Prototypes available

- **dwarf-P-cloudMicrophysics-IFSScheme**: The original cloud scheme
  from IFS that is naturally suited to host-type machines and
  optimized on the Cray system at ECMWF.
- **dwarf-cloudsc-fortran**: A cleaned up version of the CLOUDSC
  prototype that validates runs against platform and language-agnostic
  off-line reference data via HDF5 or the Serialbox package. The kernel code
  also is slightly cleaner than the original version.
- **dwarf-cloudsc-fortran-field**: A fortran version of CLOUDSC that uses Field API
  for the data structures. The intent of this version is to show how
  Field API is used in newer versions of the IFS.
- **dwarf-cloudsc-c**: Standalone C version of the kernel that has
  been generated by ECMWF tools. This relies exclusively on the Serialbox
  validation mechanism.
- **dwarf-cloudsc-gpu-kernels**: GPU-enabled version of the CLOUDSC dwarf
  that uses OpenACC and relies on the `!$acc kernels` directive to offload
  the computational kernel.
- **dwarf-cloudsc-gpu-scc**: GPU-enabled and optimized version of
  CLOUDSC that utilises the native blocked IFS memory layout via a
  "single-column coalesced" (SCC) loop layout. Here the outer NPROMA
  block loop is mapped to the OpenACC "gang" level and the kernel uses
  an inverted loop-nest where the outer horizontal loop is mapped to
  OpenACC " vector" parallelism. This variant lets the CUDA runtime
  manage temporary arrays and needs a large `PGI_ACC_CUDA_HEAPSIZE`
  (eg. `PGI_ACC_CUDA_HEAPSIZE=8GB` for 160K columns.)
- **dwarf-cloudsc-gpu-scc-hoist**: GPU-enabled and optimized version of
  CLOUDSC that also uses the SCC loop layout, but promotes the inner
  "vector" loop to the driver and declares the kernel as sequential.
  The block array arguments are fully dimensioned though, and
  multi-dimensional temporaries have been declared explicitly at the
  driver level.
- **dwarf-cloudsc-gpu-scc-k-caching**: GPU-enabled and further 
  optimized version of CLOUDSC that also uses the SCC loop layout in
  combination with loop fusion and temporary local array demotion.
- **dwarf-cloudsc-gpu-scc-cuf**: GPU-enabled and optimized version of
  CLOUDSC that uses the SCC loop layout in combination with CUDA-Fortran
  (CUF) to explicitly allocate temporary arrays in device memory and
  move parameter structures to constant memory. To enable this variant,
  a suitable CUDA installation is required and the `--with-cuda` flag
  needs to be passed at the build stage.
- **dwarf-cloudsc-gpu-scc-cuf-k-caching**: GPU-enabled and further
  optimized version of CLOUDSC that uses the SCC loop layout in
  combination with loop fusion and temporary local array demotion, implemented
  using CUDA-Fortran (CUF). To enable this variant,
  a suitable CUDA installation is required and the `--with-cuda` flag
  needs to be passed at the build stage.
- **CUDA C prototypes**: To enable these variants, a suitable 
  CUDA installation is required and the `--with-cuda` flag needs
  to be pased at the build stage.
  - **dwarf-cloudsc-cuda**: GPU-enabled, CUDA C version of CLOUDSC.
  - **dwarf-cloudsc-cuda-hoist**: GPU-enabled, optimized CUDA C version 
    of CLOUDSC including host side hoisted temporary local variables.
  - **dwarf-cloudsc-cuda-k-caching**: GPU-enabled, further optimized CUDA
    C version of CLOUDSC including loop fusion and temporary local 
    array demotion. 
  - **dwarf-cloudsc-cuda-opt**: GPU-enabled, further optimized beyond
    k-caching CUDA C version that buffers some variables and
    uses pipelined global-to-shared memory copies that are overlapped
    with compute (TMA loads). 
- **dwarf-cloudsc-gpu-scc-field**: GPU-enabled and optimized version of
  CLOUDSC that uses the SCC loop layout, and uses [FIELD API](https://github.com/ecmwf-ifs/field_api) (a Fortran library purpose-built for IFS data-structures that facilitates the
  creation and management of field objects in scientific code) to perform device offload 
  and copyback.
  The field api variant supports modern features of the FIELD API such as *field gangs* that group
  multiple fields and allocates them in one larger field, in order to reduce allocations and
  data transfers. Field gang support can be enabled at runtime by setting the environment
  variable `CLOUDSC_PACKED_STORAGE=ON`. If CUDA is available, then the field api variant also supports
  the use of allocating fields in pinned memory. This is enabled by setting the
  environemnt variable `CLOUDSC_FIELD_API_PINNED=ON` and will speed up data transfers between host and device.
  To enable this variant, a suitable CUDA installation is required and the
  `--with-cuda` flag needs to be passed at the build stage. This variant lets the CUDA runtime 
  manage temporary arrays and needs a large `NV_ACC_CUDA_HEAPSIZE` (eg. `NV_ACC_CUDA_HEAPSIZE=8GB` for 160K columns.).
  It is possible to disable Field API registering fields in the OpenACC data map, by passing the
  `--without-mapped-fields` flag at build stage.
- **cloudsc-pyiface.py**: a combination of the cloudsc/cloudsc-driver routines
  of cloudsc-fortran with the uppermost `dwarf` program replaced with a
  corresponding Python script capable of HDF5 data load and 
  verification of computation results. The computation is realized by the
  Fortran subprogram, mimicking cloudsc-fortran and equipped with only 
  minor modifications (i.e. derived types/global paramters handling).
  Turned off by default, activate at the build stage with 
  `--cloudsc-fortran-pyiface=ON`.
- **dwarf-cloudsc-fortran-atlas**: A version of **dwarf-cloudsc-fortran** which uses the [Atlas library](https://github.com/ecmwf/atlas) 
  and its Field and FieldSet data stuctures. There are two storage settings for variables. If the environment variable
  CLOUDSC_ATLAS_MULTIFIELD is "0", "OFF", or "FALSE", the variables are managed as atlas::FieldSet, which is an array of atlas::Fields. For other values of CLOUDSC_ATLAS_MULTIFIELD, a batching of variables is used as (BLK_IDX, LEV, VAR_ID, BLK_ID).
  

## Download and Installation

The code is written in Fortran 2003 and it has been tested using the various compilers, including:

    GCC 7.3, 9.3, 11.2
    Cray 8.7.7
    NVHPC 20.9, 22.1
    Intel (classic)

This application does not need MPI nor BLAS libraries for performance. Just a compiler that understands
OpenMP directives. Fortran must be at least level F2003.

Inside the dwarf directory you can find some example of outputs inside the example-outputs/ directory.

In addition, to run the dwarf it is necessary to use an input file that can be found inside the config-files/
directory winthin the dwarf folder.

The preferred method to install the CLOUDSC dwarf uses the bundle
definition shipped in the main repository. For this please
install the bundle via:

```sh
./cloudsc-bundle create  # Checks out dependency packages
./cloudsc-bundle build [--build-type=debug|bit|release] [--arch=./arch/ecmwf/machine/compiler/version]
```

The individual prototype variants of the dwarf are managed as ECBuild features
and can be enable or disabled via `--cloudsc-<feature>=[ON|OFF]` arguments to
`cloudsc-bundle build`.

The use of the `boost` library or module is required by the Serialbox
utility package for filesystem utilities. If `boost` is not available
on a given system, Serialbox's internal "experimental filesystem" can
be used via the `--serialbox-experimental=ON` argument, although this
has proven difficult with certain compiler toolchains.

### GPU versions of CLOUDSC

The GPU-enabled versions of the dwarf are by default disabled. To
enable them use the `--with-gpu` flag. For example to build on the ECMWF's ATOS
A100 nodes:

```sh
./cloudsc-bundle create  # Checks out dependency packages
./cloudsc-bundle build --clean --with-gpu --arch=./arch/ecmwf/hpc2020/nvhpc/22.1
```

### MPI-enabled versions of CLOUDSC

Optionally, dwarf-cloudsc-fortran and the GPU versions can be built with
MPI support by providing the `--with-mpi` flag. For example on ATOS:

```sh
./cloudsc-bundle create
./cloudsc-bundle build --clean --with-mpi --with-gpu --arch=./arch/ecmwf/hpc2020/nvhpc/22.1
```

Running with MPI parallelization distributes the columns of the working set
among all ranks. The specified number of OpenMP threads is then spawned on
each rank. Results are gathered from all ranks and reported for the global
working set. Performance numbers are also gathered and reported per thread,
per rank and total.

**Important:** If the total size of the working set (2nd argument, see
"[Running and testing](#running-and-testing)") **exceeds** the number of
columns in the input file (the input data in the repository consists of just
100 columns), every rank derives its working set by replicating the columns in
the input file, starting with the first column in the file. This means, all
ranks effectively work on the same data set.
If the total size of the working set is **less than or equal** to the number of
columns in the input file, these are truly distributed and every rank ends up
with a different working set.

When running with multiple GPUs each rank needs to be assigned a different
device. This can be achieved using the `CUDA_VISIBLE_DEVICES` environment
variable:

```sh
mpirun -np 2 bash -c "CUDA_VISIBLE_DEVICES=\${OMPI_COMM_WORLD_RANK} bin/dwarf-cloudsc-gpu-scc-stack 1 163840 128"
```

### Choosing between HDF5 and Serialbox input file format

The default build configuration relies on HDF5 input and reference data for
dwarf-cloudsc-fortran as well as GPU and Loki versions. The original
dwarf-P-cloudMicrophysics-IFSScheme always uses raw Fortran binary format.

**Please note:** The HDF55 installation needs to have the f03 interfaces installed (default with HDF5 1.10+).

As an alternative to HDF5, the [Serialbox](https://github.com/GridTools/serialbox)
library can be used to load input and reference data. This, however, requires
certain boost libraries or its own internal experimental filesystem, both of
which proved difficult on certain compiler toolchains or more exotic hardware
architectures.

The original input is provided as raw Fortran binary in prototype1, but
input and reference data can be regenerated from this variant by running

```sh
CLOUDSC_WRITE_INPUT=1 ./bin/dwarf-P-cloudMicrophysics-IFSScheme 1 100 100
CLOUDSC_WRITE_REFERENCE=1 ./bin/dwarf-P-cloudMicrophysics-IFSScheme 1 100 100
```

Note that this is only available via Serialbox at the moment. Updates to HDF5
input or reference data have to be done via manual conversion. A small
Python script for this with usage instructions can be found in the
[serialbox2hdf5](serialbox2hdf5/README.md) directory.

### Building on ECMWF's Atos BullSequana XH2000

To build on ECMWF's Atos BullSequana XH2000 supercomputer, run the following commands:

```sh
./cloudsc-bundle create
./cloudsc-bundle build --arch arch/ecmwf/hpc2020/compiler/version [--single-precision] [--with-mpi]
```

Currently available `compiler/version` selections are:

* `gnu/9.3.0` and `gnu/11.2.0`
* `intel/2021.4.0`
* `nvhpc/22.1` (use with `--with-gpu` on AC's GPU partition)

### A64FX version of CLOUDSC

Preliminary results for CLOUDSC have been generated for A64FX CPUs on
Isambard. A set of arch and toolchain files and detailed installation
and run instructions are provided
[here](https://confluence.ecmwf.int/display/~nabr/3rd+Isambard+Hackathon).

### SYCL version of CLOUDSC

A preliminary SYCL code variant has been added and tested with a custom
DPCPP install on ECMWF's AC partition. To build this, please use the
SYCL-specific environment setups:

```
./cloudsc-bundle build --clean --build-dir=build-sycl --with-gpu --with-sycl --with-serialbox --arch=arch/ecmwf/hpc2020/intel-sycl/2021.4.0

# Then run with
cd build-sycl && . env.sh
./bin/dwarf-cloudsc-scc-sycl 1 240000 128
./bin/dwarf-cloudsc-scc-hoist-sycl 1 240000 128
./bin/dwarf-cloudsc-scc-k-caching-sycl 1 240000 128
```

## Running and testing

The different prototype variants of the dwarf create different binaries that
all behave similarly. The basic three arguments define (in this order):

- Number of OpenMP threads
  - 1 : single thread mode, skip multithread MPI init; default value;
  - 2 or higher : force OpenMP thread count, enables multithread MPI;
  - 0 or negative : read OMP_NUM_THREADS variable if present or defaults to CPU count (`omp_get_max_threads()`);
- Size of overall working set in columns
- Block size (NPROMA) in columns

An example:

```sh
cd build
./bin/dwarf-P-cloudMicrophysics-IFSScheme 4 16384 32  # The original
./bin/dwarf-cloudsc-fortran 4 16384 32   # The cleaned-up Fortran
./bin/dwarf-cloudsc-c 4 16384 32   # The standalone C version
```

### Running on ECMWF's Atos BullSequana XH2000

On the Atos system, a high-watermark run on a single socket can be performed as follows:

```sh
export OMP_NUM_THREADS=64
OMP_PLACES="{$(seq -s '},{' 0 $(($OMP_NUM_THREADS-1)) )}" srun -q np --ntasks=1 --hint=nomultithread --cpus-per-task=$OMP_NUM_THREADS ./bin/dwarf-cloudsc-fortran $OMP_NUM_THREADS 163840 32
```

For a double-precision build with the GNU 11.2.0 compiler, performance of
~73 GF/s is achieved.

To run the GPU variant on AC, which includes some GPU nodes, allocate
an interactive session on a GPU node and run the binary as usual:

```sh
srun -N1 -q ng -p gpu --gres=gpu:4 --mem 200G --pty /bin/bash
bin/dwarf-cloudsc-gpu-scc-hoist 1 262144 128
```

For a double-precision build with NVHPC 22.1, performance of ~340 GF/s
on a single GPU is achieved.

A multi-GPU run requires MPI (build with `--with-mpi`) with a dedicated MPI
task for each GPU and (at the moment) manually assigning CUDA devices to each
rank, as Slurm is not yet fully configured for the GPU partition.

To use four GPUs on one node, allocate the relevant resources
```sh
salloc -N 1 --tasks-per-node 4 -q ng -p gpu --gres=gpu:4 --mem 200G
```

and then run the binary like this:

```sh
srun bash -c "CUDA_VISIBLE_DEVICES=\$SLURM_LOCALID bin/dwarf-cloudsc-gpu-scc-hoist 1 \$((\$SLURM_NPROCS*262144)) 128"
```

In principle, the same should work for multi-node execution (`-N 2`, `-N 4` etc.) once interconnect issues are resolved.

### GPU runs: Timing device kernels and data transfers

For GPU-enabled runs two internal timer results are reported:

* The isolated compute time of the main compute kernel on device (where `#BLKS == 1`)
* The overall time of the execution loop including data offload and copyback

It is important to note that due to the nature of the kernel, data
transfer overheads will dominate timings, and that most supported GPU
variants aim to optimise compute kernel timings only. However, a
dedicated variant `dwarf-cloudsc-gpu-scc-field` has been added to
explore host-side memory pinning, which improves data transfer times
and alternative data layout strategies. By default, pinned memory is turned off
but can be turned on by setting the environment variable `CLOUDSC_FIELD_API_PINNED=ON`.
This will allocate each array variable individually in pinned memory. A runtime flag
`CLOUDSC_PACKED_STORAGE=ON` can be used to enable "packed" storage,
where multiple arrays are stored in a single base allocation, eg.

```sh
NV_ACC_CUDA_HEAPSIZE=8G CLOUDSC_PACKED_STORAGE=ON ./bin/dwarf-cloudsc-gpu-scc-field 1 80000 128
```

## Loki transformations for CLOUDSC

[Loki](https://github.com/ecmwf-ifs/loki) is an in-house developed
source-to-source translation tool that allows us to create bespoke
transformations for the IFS to target and experiment with emerging HPC
architectures and programming models. We use the CLOUDSC dwarf as a demonstrator
for targeted transformation capabilities of physics and grid point computations
kernels, including conversion to C and GPU.

The following build flags enable the demonstrator build targets on the
ECMWF Atos HPC facility's GPU partition:

```sh
./cloudsc-bundle build --clean [--with-gpu] --with-loki --loki-frontend=fp --arch=./arch/ecmwf/hpc2020/nvhpc/22.1
```

The following Loki modes are included in the dwarf, each with a bespoke demonstrator build:

- **cloudsc-loki-idem**: "Idempotence" mode that performs a full
  parse-unparse cycle of the kernel and performs various housekeeping
  transformations, including the driver-level source injection
  mechanism currently facilitated by Loki.
- **cloudsc-loki-sca**: Pure single-column mode that strips all horizontal
  vector loops from the kernel and introduces an outer "column-loop"
  at the driver level.
- **cloudsc-loki-c**: A prototype C transpilation pipeline that converts
  the kernel to C and calls it via iso_c_bindings interfaces from the
  driver.

## Python-driven CLOUDSC variants
The following partly or fully Python-based CLOUDSC are available:
- **cloudsc-python**: GT4PY based Python-only implementation. Refer to `src/cloudsc_python`
  for information on how to bootstrap/execute this variant
- **cloudsc-pyiface**: Fortran-based CLOUDSC variant driven by the Python script.
  Activate with: 
```sh
./cloudsc-bundle build --clean --cloudsc-fortran-pyiface=ON
```
These variants are disabled by default. Refer to README.md in corresponding subdirectories
for further information.
  
### A note on frontends

Loki currently supports three frontends to parse the Fortran source code:

- [FParser](https://github.com/stfc/fparser) (`loki-frontend=fp`):
  The preferred default; developed by STFC for PsyClone.
- [OMNI](https://github.com/omni-compiler/omni-compiler) frontend (`loki-frontend=omni`):

For completeness, all three frontends are tested in our CI, which
means we require the `.xmod` module description files for utility
routines in `src/common` for processing the CLOUDSC source files with
the OMNI frontend. These are stored in the source under
`src/cloudsc_loki/xmod`.

## Benchmarking

To automate parameter space sweeps and ease testing across various platforms, a
[JUBE](https://www.fz-juelich.de/jsc/jube) benchmark definition is included in
the directory `benchmark`. See the included [README](benchmark/README.md) for
further details and usage instructions.
