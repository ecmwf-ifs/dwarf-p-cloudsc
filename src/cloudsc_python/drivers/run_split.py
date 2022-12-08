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
from typing import Optional

from cloudsc4py.physics.cloudsc_split import Cloudsc

from config import default_python_config, default_io_config
from run import core


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
    host_alias: Optional[str],
    output_csv_file: Optional[str],
    output_csv_file_stencils: Optional[str],
) -> None:
    """
    Driver for the GT4Py-based implementation of CLOUDSC.

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
    core(config, io_config, cloudsc_cls=Cloudsc)

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


if __name__ == "__main__":
    main()
