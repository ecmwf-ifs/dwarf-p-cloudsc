# JUBE benchmarking environment

JUBE is a benchmarking environment that provides a script-based framework to
create benchmark sets, run them and evaluate the results. It is developed at
Forschungszentrum Juelich, Germany.

Further information: https://www.fz-juelich.de/jsc/jube

Documentation: https://apps.fz-juelich.de/jsc/jube/jube2/docu/

## Installation

```bash
# Python 3 module loaded
python3 -m venv venv
source venv/bin/activate
pip install -r requirements.txt
```

## Running a benchmark

```bash
# Virtual environment loaded
jube run <benchmark>.yml
```

## Analysing and displaying results

```bash
# Virtual environment loaded
jube analyse <benchmark run directory> --id <benchmark id>
jube result <benchmark run directory> --id <benchmark id> | less -S
```

or, both in one:

```bash
# Virtual environment loaded
jube result -a <benchmark run directory> --id <benchmark id> | less -S
```

Skipping the benchmark id is equivalent to using the latest benchmark run.

## Useful commands

Update benchmark results after modifying patterns or result table, without
re-running the benchmark:

```bash
# Virtual environment loaded
jube result -a -u <benchmark>.yml <benchmark run directory> --id <benchmark id> | less -S
```

Run without actually executing benchmark steps (useful to validate the YAML
files and see the parameter space expansion):

```bash
# Virtual environment loaded
jube --debug -v run <benchmark>.yml
```
