# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.
# (C) Copyright 2022- ETH Zurich.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import click
import os
import subprocess
from typing import Optional

from config import FortranConfig, IOConfig, default_fortran_config, default_io_config
from utils import print_performance, to_csv


def core(config: FortranConfig, io_config: IOConfig) -> None:
    executable = os.path.join(
        os.path.dirname(__file__), config.build_dir, f"bin/dwarf-cloudsc-{config.variant}"
    )
    if not os.path.exists(executable):
        raise RuntimeError(f"The executable `{executable}` does not exist.")

    # warm-up cache
    _ = subprocess.run(
        [
            executable,
            str(config.num_threads),
            str(config.num_cols),
            str(min(config.num_cols, config.nproma)),
        ],
        capture_output=True,
    )

    # run and profile
    runtimes = []
    for _ in range(config.num_runs):
        out = subprocess.run(
            [
                executable,
                str(config.num_threads),
                str(config.num_cols),
                str(min(config.num_cols, config.nproma)),
            ],
            capture_output=True,
        )
        if "gpu" in config.variant:
            x = out.stderr.decode("utf-8").split("\n")[2]
            y = x.split(" ")
            z = [c for c in y if c != ""]
            runtimes.append(float(z[-4]))
        else:
            x = out.stderr.decode("utf-8").split("\n")[-2]
            y = x.split(" ")
            z = [c for c in y if c != ""]
            runtimes.append(float(z[-4]))

    runtime_mean, runtime_stddev = print_performance(runtimes)

    if io_config.output_csv_file is not None:
        to_csv(
            io_config.output_csv_file,
            io_config.host_name,
            config.variant,
            config.num_cols,
            config.num_runs,
            runtime_mean,
            runtime_stddev,
        )


@click.command()
@click.option(
    "--build-dir",
    type=str,
    default="fortran",
    help="Path to the build directory of the FORTRAN dwarf.",
)
@click.option(
    "--variant",
    type=str,
    default="fortran",
    help="Code variant."
    "\n\nOptions: fortran, gpu-scc, gpu-scc-hoist, gpu-omp-scc-hoist."
    "\n\nDefault: fortran.",
)
@click.option(
    "--nproma",
    type=int,
    default=32,
    help="Block size.\n\nRecommended values: 32 on CPUs, 128 on GPUs.\n\nDefault: 32.",
)
@click.option("--num-cols", type=int, default=1, help="Number of domain columns.\n\nDefault: 1.")
@click.option(
    "--num-runs",
    type=int,
    default=1,
    help="Number of executions.\n\nDefault: 1.",
)
@click.option(
    "--num-threads",
    type=int,
    default=1,
    help="Number of threads."
    "\n\nRecommended values: 24 on Piz Daint's CPUs, 128 on MLux's CPUs, 1 on GPUs."
    "\n\nDefault: 1.",
)
@click.option("--host-alias", type=str, default=None, help="Name of the host machine (optional).")
@click.option(
    "--output-csv-file",
    type=str,
    default=None,
    help="Path to the CSV file where writing performance counters (optional).",
)
def main(
    build_dir: str,
    variant: str,
    nproma: int,
    num_cols: int,
    num_runs: int,
    num_threads: int,
    host_alias: Optional[str],
    output_csv_file: Optional[str],
) -> None:
    """Driver for the FORTRAN implementation of CLOUDSC."""
    config = (
        default_fortran_config.with_build_dir(build_dir)
        .with_variant(variant)
        .with_nproma(nproma)
        .with_num_cols(num_cols)
        .with_num_runs(num_runs)
        .with_num_threads(num_threads)
    )
    io_config = default_io_config.with_output_csv_file(output_csv_file).with_host_name(host_alias)
    core(config, io_config)


if __name__ == "__main__":
    main()
