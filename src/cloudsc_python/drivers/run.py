# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.
# (C) Copyright 2022- ETH Zurich.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import click
import csv
import datetime
import numpy as np
import os
from typing import Optional, Type
import subprocess

from cloudsc4py.framework.grid import ComputationalGrid
from cloudsc4py.framework.config import DataTypes, GT4PyConfig, PythonConfig, IOConfig, FortranConfig
from cloudsc4py.physics.cloudsc import Cloudsc
from cloudsc4py.physics.cloudsc_split import CloudscSplit
from cloudsc4py.initialization.reference import get_reference_tendencies, get_reference_diagnostics
from cloudsc4py.initialization.state import get_state
from cloudsc4py.utils.iox import HDF5Reader
from cloudsc4py.utils.timing import timing
from cloudsc4py.utils.validation import validate


# Default configurations for the various run modes

default_io_config = IOConfig(output_file=None, host_name=None)

config_files_dir = os.path.normpath(os.path.join(os.path.dirname(__file__), "../../../config-files"))
default_python_config = PythonConfig(
    num_cols=1,
    enable_validation=True,
    input_file=os.path.join(config_files_dir, "input.h5"),
    reference_file=os.path.join(config_files_dir, "reference.h5"),
    num_runs=15,
    data_types=DataTypes(bool=bool, float=np.float64, int=int),
    gt4py_config=GT4PyConfig(backend="numpy", rebuild=False, validate_args=True, verbose=True),
    sympl_enable_checks=True,
)


default_fortran_config = FortranConfig(
    build_dir=".", variant="fortran", nproma=32, num_cols=1, num_runs=1, num_threads=1
)


# Some utility function for performance benchmarking

def to_csv(
    output_file: str,
    host_name: str,
    variant: str,
    num_cols: int,
    num_runs: int,
    runtime_mean: float,
    runtime_stddev: float,
) -> None:
    """Write mean and standard deviation of measured runtimes to a CSV file."""
    if not os.path.exists(output_file):
        with open(output_file, "w") as csv_file:
            writer = csv.writer(csv_file, delimiter=",")
            writer.writerow(("date", "host", "variant", "num_cols", "num_runs", "mean", "stddev"))
    with open(output_file, "a") as csv_file:
        writer = csv.writer(csv_file, delimiter=",")
        writer.writerow(
            (
                datetime.date.today().strftime("%Y%m%d"),
                host_name,
                variant,
                num_cols,
                num_runs,
                runtime_mean,
                runtime_stddev,
            )
        )


def print_performance(runtimes):
    """Print means and standard deviation of measure runtimes to screen."""
    n = len(runtimes)
    mean = sum(runtimes) / n
    stddev = (sum((runtime - mean) ** 2 for runtime in runtimes) / (n - 1 if n > 1 else n)) ** 0.5
    print(f"Performance: Average runtime over {n} runs: {mean:.3f} \u00B1 {stddev:.3f} ms.")
    return mean, stddev


# The core benchmark implementations in GT4Py and for external Fortran runners

def core_cloudsc(config: PythonConfig, io_config: IOConfig, cloudsc_cls: Type) -> None:
    hdf5_reader = HDF5Reader(config.input_file, config.data_types)

    nx = config.num_cols or hdf5_reader.get_nlon()
    nz = hdf5_reader.get_nlev()
    computational_grid = ComputationalGrid(nx, 1, nz)

    state = get_state(computational_grid, hdf5_reader, gt4py_config=config.gt4py_config)
    dt = hdf5_reader.get_timestep()

    yoecldp_paramaters = hdf5_reader.get_yoecldp_parameters()
    yoethf_parameters = hdf5_reader.get_yoethf_parameters()
    yomcst_parameters = hdf5_reader.get_yomcst_parameters()
    yrecldp_parameters = hdf5_reader.get_yrecldp_parameters()

    cloudsc = cloudsc_cls(
        computational_grid,
        yoecldp_paramaters,
        yoethf_parameters,
        yomcst_parameters,
        yrecldp_parameters,
        enable_checks=config.sympl_enable_checks,
        gt4py_config=config.gt4py_config,
    )
    tends, diags = cloudsc(state, dt)

    runtimes = []
    for i in range(config.num_runs):
        with timing(f"run_{i}") as timer:
            cloudsc(state, dt, out_tendencies=tends, out_diagnostics=diags)
        runtimes.append(timer.get_time(f"run_{i}") * 1000)

    runtime_mean, runtime_stddev = print_performance(runtimes)

    if io_config.output_csv_file is not None:
        to_csv(
            io_config.output_csv_file,
            io_config.host_name,
            config.gt4py_config.backend,
            nx,
            config.num_runs,
            runtime_mean,
            runtime_stddev,
        )

    if config.enable_validation:
        hdf5_reader_ref = HDF5Reader(config.reference_file, config.data_types)
        tends_ref = get_reference_tendencies(
            computational_grid, hdf5_reader_ref, gt4py_config=config.gt4py_config
        )
        diags_ref = get_reference_diagnostics(
            computational_grid, hdf5_reader_ref, gt4py_config=config.gt4py_config
        )

        tends_fail = validate(tends, tends_ref)
        if len(tends_fail) == 0:
            print("Results: All tendencies have been successfully validated. HOORAY!")
        else:
            print(
                f"Results: Validation failed for {len(tends_fail)}/{len(tends_ref) - 1} "
                f"tendencies: {', '.join(tends_fail)}."
            )

        diags_fail = validate(diags, diags_ref)
        if len(diags_fail) == 0:
            print("Results: All diagnostics have been successfully validated. HOORAY!")
        else:
            print(
                f"Results: Validation failed for {len(diags_fail)}/{len(diags_ref) - 1} "
                f"diagnostics: {', '.join(diags_fail)}."
            )


def core_fortran(config: FortranConfig, io_config: IOConfig) -> None:
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


# Click group to present various CLI modes under a single runner script

@click.group()
def main():
    """
    Driver for the GT4Py-based implementation of CLOUDSC.
    """
    pass


@main.command()
@click.option(
    "--backend", type=str, default=None, help="GT4Py backend."
    "\n\nOptions: numpy, gt:cpu_kfirst, gt:cpu_ifirst, gt:gpu, cuda, dace:cpu, dace:gpu."
    "\n\nDefault: numpy.",
)
@click.option(
    "--enable-checks/--disable-checks", is_flag=True, type=bool, default=False,
    help="Enable/disable sanity checks performed by Sympl and GT4Py.\n\nDefault: enabled.",
)
@click.option(
    "--enable-validation/--disable-validation", is_flag=True, type=bool, default=True,
    help="Enable/disable data validation.\n\nDefault: enabled.",
)
@click.option(
    "--num-cols", type=int, default=None,
    help="Number of domain columns.\n\nDefault: 1."
)
@click.option(
    "--num-runs", type=int, default=1,
    help="Number of executions.\n\nDefault: 1.",
)
@click.option(
    "--host-alias", type=str, default=None,
    help="Name of the host machine (optional)."
)
@click.option(
    "--output-csv-file", type=str, default=None,
    help="Path to the CSV file where writing performance counters (optional).",
)
@click.option(
    "--output-csv-file-stencils", type=str, default=None,
    help="Path to the CSV file where writing performance counters for each stencil (optional).",
)
def single(
    backend: Optional[str],
    enable_checks: bool,
    enable_validation: bool,
    num_cols: Optional[int],
    num_runs: Optional[int],
    host_alias: Optional[str],
    output_csv_file: Optional[str],
    output_csv_file_stencils: Optional[str],
) -> None:
    """
    Computations are carried out in a single stencil.
    """
    config = (
        default_python_config.with_backend(backend)
        .with_checks(enable_checks)
        .with_validation(enable_validation)
        .with_num_cols(num_cols)
        .with_num_runs(num_runs)
    )
    io_config = default_io_config.with_output_csv_file(output_csv_file).with_host_name(host_alias)
    core_cloudsc(config, io_config, cloudsc_cls=Cloudsc)

    if output_csv_file_stencils is not None:
        call_time = None
        for key, value in config.gt4py_config.exec_info.items():
            if "cloudsc" in key:
                call_time = value["total_call_time"] * 1000 / config.num_runs

        if not os.path.exists(output_csv_file_stencils):
            with open(output_csv_file_stencils, "w") as f:
                writer = csv.writer(f, delimiter=",")
                writer.writerow(("date", "host", "backend", "num_cols", "num_runs", "cloudsc"))
        with open(output_csv_file_stencils, "a") as f:
            writer = csv.writer(f, delimiter=",")
            writer.writerow(
                (
                    datetime.date.today().strftime("%Y%m%d"),
                    io_config.host_name,
                    config.gt4py_config.backend,
                    config.num_cols,
                    config.num_runs,
                    call_time,
                )
            )


@main.command()
@click.option(
    "--backend",  type=str, default=None, help="GT4Py backend."
    "\n\nOptions: numpy, gt:cpu_kfirst, gt:cpu_ifirst, gt:gpu, cuda, dace:cpu, dace:gpu."
    "\n\nDefault: numpy.",
)
@click.option(
    "--enable-checks/--disable-checks", is_flag=True, type=bool, default=False,
    help="Enable/disable sanity checks performed by Sympl and GT4Py.\n\nDefault: enabled.",
)
@click.option(
    "--enable-validation/--disable-validation", is_flag=True, type=bool, default=True,
    help="Enable/disable data validation.\n\nDefault: enabled.",
)
@click.option(
    "--num-cols", type=int, default=None,
    help="Number of domain columns.\n\nDefault: 1."
)
@click.option(
    "--num-runs", type=int, default=1,
    help="Number of executions.\n\nDefault: 1.",
)
@click.option(
    "--host-alias", type=str, default=None,
    help="Name of the host machine (optional)."
)
@click.option(
    "--output-csv-file", type=str, default=None,
    help="Path to the CSV file where writing performance counters (optional).",
)
@click.option(
    "--output-csv-file-stencils", type=str, default=None,
    help="Path to the CSV file where writing performance counters for each stencil (optional).",
)
def split(
    backend: Optional[str],
    enable_checks: bool,
    enable_validation: bool,
    num_cols: Optional[int],
    num_runs: Optional[int],
    host_alias: Optional[str],
    output_csv_file: Optional[str],
    output_csv_file_stencils: Optional[str],
) -> None:
    """
    Computations are split into two stencils.
    """
    config = (
        default_python_config.with_backend(backend)
        .with_checks(enable_checks)
        .with_validation(enable_validation)
        .with_num_cols(num_cols)
        .with_num_runs(num_runs)
    )
    io_config = default_io_config.with_output_csv_file(output_csv_file).with_host_name(host_alias)
    core_cloudsc(config, io_config, cloudsc_cls=CloudscSplit)

    if output_csv_file_stencils is not None:
        cloudsc_tendencies_call_time = None
        cloudsc_fluxes_call_time = None
        for key, value in config.gt4py_config.exec_info.items():
            if "tendencies" in key:
                cloudsc_tendencies_call_time = value["total_call_time"] * 1000 / config.num_runs
            elif "fluxes" in key:
                cloudsc_fluxes_call_time = value["total_call_time"] * 1000 / config.num_runs

        if not os.path.exists(output_csv_file_stencils):
            with open(output_csv_file_stencils, "w") as f:
                writer = csv.writer(f, delimiter=",")
                writer.writerow(
                    (
                        "date",
                        "host",
                        "backend",
                        "num_cols",
                        "num_runs",
                        "cloudsc_tendencies",
                        "cloudsc_fluxes",
                    )
                )
        with open(output_csv_file_stencils, "a") as f:
            writer = csv.writer(f, delimiter=",")
            writer.writerow(
                (
                    datetime.date.today().strftime("%Y%m%d"),
                    io_config.host_name,
                    config.gt4py_config.backend,
                    config.num_cols,
                    config.num_runs,
                    cloudsc_tendencies_call_time,
                    cloudsc_fluxes_call_time,
                )
            )


@main.command()
@click.option(
    "--build-dir", type=str, default="fortran",
    help="Path to the build directory of the FORTRAN dwarf.",
)
@click.option(
    "--variant", type=str, default="fortran", help="Code variant."
    "\n\nOptions: fortran, gpu-scc, gpu-scc-hoist, gpu-omp-scc-hoist."
    "\n\nDefault: fortran.",
)
@click.option(
    "--nproma", type=int, default=32,
    help="Block size.\n\nRecommended values: 32 on CPUs, 128 on GPUs.\n\nDefault: 32.",
)
@click.option(
    "--num-cols", type=int, default=1,
    help="Number of domain columns.\n\nDefault: 1."
)
@click.option(
    "--num-runs", type=int, default=1,
    help="Number of executions.\n\nDefault: 1.",
)
@click.option(
    "--num-threads", type=int, default=1, help="Number of threads."
    "\n\nRecommended values: 24 on Piz Daint's CPUs, 128 on MLux's CPUs, 1 on GPUs."
    "\n\nDefault: 1.",
)
@click.option(
    "--host-alias", type=str, default=None,
    help="Name of the host machine (optional)."
)
@click.option(
    "--output-csv-file", type=str, default=None,
    help="Path to the CSV file where writing performance counters (optional).",
)
def fortran(
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
    core_fortran(config, io_config)


if __name__ == "__main__":
    main()
