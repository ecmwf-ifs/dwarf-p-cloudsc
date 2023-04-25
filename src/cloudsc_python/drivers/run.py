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
import os
from typing import Optional, Type

from cloudsc4py.framework.grid import ComputationalGrid
from cloudsc4py.physics.cloudsc import Cloudsc
from cloudsc4py.initialization.reference import get_reference_tendencies, get_reference_diagnostics
from cloudsc4py.initialization.state import get_state
from cloudsc4py.utils.iox import HDF5Reader
from cloudsc4py.utils.timing import timing
from cloudsc4py.utils.validation import validate

from config import PythonConfig, IOConfig, default_python_config, default_io_config
from utils import print_performance, to_csv


def core(config: PythonConfig, io_config: IOConfig, cloudsc_cls: Type) -> None:
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
        runtimes.append(timer.get_time(f"run_{i}"))

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


@click.command()
@click.option(
    "--backend",
    type=str,
    default=None,
    help="GT4Py backend."
    "\n\nOptions: numpy, gt:cpu_kfirst, gt:cpu_ifirst, gt:gpu, cuda, dace:cpu, dace:gpu."
    "\n\nDefault: numpy.",
)
@click.option(
    "--enable-checks/--disable-checks",
    is_flag=True,
    type=bool,
    default=False,
    help="Enable/disable sanity checks performed by Sympl and GT4Py.\n\nDefault: enabled.",
)
@click.option(
    "--enable-validation/--disable-validation",
    is_flag=True,
    type=bool,
    default=True,
    help="Enable/disable data validation.\n\nDefault: enabled.",
)
@click.option("--num-cols", type=int, default=None, help="Number of domain columns.\n\nDefault: 1.")
@click.option(
    "--num-runs",
    type=int,
    default=1,
    help="Number of executions.\n\nDefault: 1.",
)
@click.option(
    "--precision",
    type=str,
    default="double",
    help="Select either `double` (default) or `single` precision.",
)
@click.option("--host-alias", type=str, default=None, help="Name of the host machine (optional).")
@click.option(
    "--output-csv-file",
    type=str,
    default=None,
    help="Path to the CSV file where writing performance counters (optional).",
)
@click.option(
    "--output-csv-file-stencils",
    type=str,
    default=None,
    help="Path to the CSV file where writing performance counters for each stencil (optional).",
)
def main(
    backend: Optional[str],
    enable_checks: bool,
    enable_validation: bool,
    num_cols: Optional[int],
    num_runs: Optional[int],
    precision: str,
    host_alias: Optional[str],
    output_csv_file: Optional[str],
    output_csv_file_stencils: Optional[str],
) -> None:
    """
    Driver for the GT4Py-based implementation of CLOUDSC.

    Computations are carried out in a single stencil.
    """
    config = (
        default_python_config.with_backend(backend)
        .with_checks(enable_checks)
        .with_validation(enable_validation)
        .with_num_cols(num_cols)
        .with_num_runs(num_runs)
        .with_precision(precision)
    )
    io_config = default_io_config.with_output_csv_file(output_csv_file).with_host_name(host_alias)
    core(config, io_config, cloudsc_cls=Cloudsc)

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


if __name__ == "__main__":
    main()
