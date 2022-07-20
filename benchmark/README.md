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
`arch/<site>/<platform>/<toolchain>/<version>/include_arch.yml` but may also
include any of the other include files.

```bash
# Create a virtual environment and install JUBE
python3 -m venv venv
venv/bin/pip install -r requirements.txt

# Execute the benchmark with the correct architecture file
venv/bin/jube run cloudsc.yml --include arch/<site>/<platform>/<toolchain>/<version> [--only-bench=<cpu|gpu>]

# Analyse output and create results table
venv/bin/jube result -a rundir_<cpu|gpu> --id=<benchmark id> | less -S
```
