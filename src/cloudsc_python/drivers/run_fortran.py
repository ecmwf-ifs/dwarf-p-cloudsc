# -*- coding: utf-8 -*-
import click
import os
import subprocess
from typing import Optional

from config import FortranConfig, IOConfig, default_fortran_config, default_io_config
from utils import print_performance, to_csv


def core(config: FortranConfig, io_config: IOConfig) -> None:
    executable = os.path.join(
        os.path.dirname(__file__), config.build_dir, f"dwarf-cloudsc-{config.mode}"
    )
    if not os.path.exists(executable):
        raise RuntimeError(f"The executable `{executable}` does not exist.")

    # warm-up cache
    _ = subprocess.run(
        [
            executable,
            str(config.num_threads),
            str(config.nx),
            str(min(config.nx, 32)),
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
                str(config.nx),
                str(min(config.nx, 32)),
            ],
            capture_output=True,
        )
        x = out.stderr.decode("utf-8").split("\n")[-2]
        y = x.split(" ")
        z = [c for c in y if c != ""]
        runtimes.append(float(z[-4]))

    runtime_mean, runtime_stddev = print_performance(runtimes)

    if io_config.output_file is not None:
        to_csv(
            io_config.output_file,
            io_config.host_name,
            config.mode,
            config.nx,
            config.num_runs,
            runtime_mean,
            runtime_stddev,
        )


@click.command()
@click.option("--build-dir", type=str, default="fortran")
@click.option("--mode", type=str, default="fortran")
@click.option("--num-runs", type=int, default=1)
@click.option("--num-threads", type=int, default=1)
@click.option("--nx", type=int, default=1)
@click.option("--output-file", type=str, default=None)
@click.option("--host-alias", type=str, default=None)
def main(
    mode: str,
    num_runs: int,
    num_threads: int,
    nx: int,
    output_file: Optional[str],
    host_alias: Optional[str],
) -> None:
    config = (
        default_fortran_config.with_mode(mode)
        .with_num_runs(num_runs)
        .with_num_threads(num_threads)
        .with_nx(nx)
    )
    io_config = default_io_config.with_output_file(output_file).with_host_name(host_alias)
    core(config, io_config)


if __name__ == "__main__":
    main()
