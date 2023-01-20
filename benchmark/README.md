# CLOUDSC benchmark setup for JUBE

This provides different benchmarks:

- cpu
- gpu

They use the JUBE benchmarking environment to automate their execution for a
range of different build and run configurations. See [JUBE.md](JUBE.md) for
an introduction.

The benchmarks are defined in the main configuration file
[`cloudsc.yml`](cloudsc.yml) and use external include files that provide the
bits and pieces to configure their execution, with bespoke overrides for each
benchmark in the main configuration file. The most relevant include files for
users are:

- [`include_arch.yml`](include/include_arch.yml):
  Hardware configuration and `arch` file to use
- [`include_parameterset.yml`](include/include_parameters.yml):
  Defines the matrix of active build options
- [`include_run.yml`](include/include_run.yml):
  Defines the matrix of execution options for each build

Probably less likely to require user editing are the following files, that
implement the mechanics of the benchmark execution:

- [`include_step.yml`](include/include_step.yml):
  Implementation of the benchmark steps
- [`include_fileset_substituteset.yml`](include/include_fileset_substituteset.yml):
  Required script file templates and required substitutions for each of them
- [`include_patternset.yml`](include/include_patternset.yml):
  Regular expressions to parse output
- [`include_analyser.yml`](include/include_analyser.yml):
  Application of regex patterns to execution stdout to collect output data
- [`include_result.yml`](include/include_result.yml):
  Compilation of result tables from analyser output data

## Usage

If it does not exist, yet, create target platform specific include files that
overwrite relevant configuration values. Typically, this means at least a
bespoke copy of `include_arch.yml` in something like
`arch/<site>/<platform>/<toolchain>/<version>/include_arch.yml` but may for certain
scenarios also require customization of the other include files.

Note that the parametersets can be initialized with the default values from the
`include` directory by providing the `init_with` option. This allows to only
overwrite values that need to be changed, e.g.

```yaml
parameterset:
  - name: arch_set
    init_with: include/include_arch.yml
    parameter:
    - {name: arch, _: arch/hpc2020/gnu/9.3.0}
```

With platform-specific overrides in place, JUBE can be installed and the benchmark
executed using the following steps:

```bash
# Create a virtual environment and install JUBE
python3 -m venv venv
venv/bin/pip install -r requirements.txt

# Execute the benchmark with the correct architecture file
venv/bin/jube run cloudsc.yml --include arch/<site>/<platform>/<toolchain>/<version> \
  [--only-bench=<cpu|gpu>] [-t <tag> [-t <tag> ...]] [-m "<description>"]

# Analyse output and create results table
venv/bin/jube result -a rundir_<cpu|gpu> --id=<benchmark id> | less -S
```

Note the following options to the `run` command:

- `--include`: This should point to the directory with the platform-specific
  include files to override parameters. This takes precedence (but does not
  replace) the default include path. Multiple include paths can be specified.
- `--only-bench`: By specifying `cpu` or `gpu`, only the relevant benchmark
  variant is being executed.
- `-t`: This allows to provide a "tag" to select certain readily available
  variations of parameters. Multiple tags can be supplied, separated by
  white space. Currently available:
  - `dp`/`sp` to enable double (the default) and single precision builds
  - `serialbox` to use Serialbox instead of HDF5 as input library
  - `mpi` to build with MPI support
  - `sweep_nproma` varies the `default_nproma` value specified for the benchmark
    by running with 1/4, 1/2, 1, 2, 4 times that value to find the optimum
- `-m`: This allows to provide a description for the benchmark execution to
  help identify a specific run later on

To view information about the performed benchmark runs, use the `info` command:

```bash
venv/bin/jube info rundir_<cpu|gpu> [--id=<benchmark id>]
```

Without `--id`, this lists all past runs and includes the description provided
via `-m`. When a benchmark id is specified, it gives a summary of that specific
benchmark run.

To postprocess the result tables, the output format can be changed to CSV by
adding `-s csv` to the `result` command.

## Typical benchmarking workflow

A typical benchmarking workflow for a new CPU platform may look like this:

1. Create a platform-specific `include_arch.yml` file
2. Install JUBE
3. Run an `NPROMA` sweep for single and double precision:
   `venv/bin/jube run cloudsc.yml --only-bench=cpu --include <path/to/dir/with/include_arch.yml> -t sp dp sweep_nproma -m "<Compiler> NPROMA sweep for platform xyz"`
4. View results to select optimum NPROMA value:
   `venv/bin/jube result rundir_cpu -a | less -S`
5. If required: Create platform-specific `include_run.yml` file that changes
   `default_nproma` value to optimum. Repeat NPROMA sweep, if necessary.
6. Run benchmark (optionally with MPI across NUMA domains), e.g. as follows:
   `venv/bin/jube run cloudsc.yml --only-bench=cpu --include <path/to/dir/with/include_arch.yml> -t sp dp mpi -m "<Compiler> MPI platform xyz"`
7. View results:
   `venv/bin/jube result rundir_cpu -a | less -S`
8. Optional: dump results to CSV:
   `venv/bin/jube result rundir_cpu -s csv > results.csv`
