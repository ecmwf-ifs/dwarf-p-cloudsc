# -*- coding: utf-8 -*-
import click
from typing import Optional

from cloudsc4py.framework.grid import ComputationalGrid
from cloudsc4py.physics.cloudsc import Cloudsc
from cloudsc4py.initialization.reference import get_reference_tendencies, get_reference_diagnostics
from cloudsc4py.initialization.state import get_state
from cloudsc4py.utils.iox import HDF5Reader
from cloudsc4py.utils.numpyx import prepare_numpy
from cloudsc4py.utils.timing import timing
from cloudsc4py.utils.validation import validate

from config import Config, default_config


def core(config: Config) -> None:
    with prepare_numpy():
        hdf5_reader = HDF5Reader(config.input_file, config.data_types)

        nx = config.nx or hdf5_reader.get_nlon()
        nz = hdf5_reader.get_nlev()
        computational_grid = ComputationalGrid(nx, 1, nz)

        state = get_state(computational_grid, hdf5_reader, gt4py_config=config.gt4py_config)
        dt = hdf5_reader.get_timestep()

        yoecldp_paramaters = hdf5_reader.get_yoecldp_parameters()
        yoethf_parameters = hdf5_reader.get_yoethf_parameters()
        yomcst_parameters = hdf5_reader.get_yomcst_parameters()
        yrecldp_parameters = hdf5_reader.get_yrecldp_parameters()

        cloudsc = Cloudsc(
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
            runtimes.append(timer.get_time(f"run_{i}", units="ms"))

        runtime_avg = sum(runtimes) / len(runtimes)
        runtime_stddev = (
            sum((runtime - runtime_avg) ** 2 for runtime in runtimes) / len(runtimes)
        ) ** 0.5
        print(
            f"{config.num_runs} runs completed in {runtime_avg:.3f} \u00B1 {runtime_stddev:.3f} ms."
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
            print("All tendencies have been successfully validated. HOORAY!")
        else:
            print(
                f"Validation failed for {len(tends_fail)}/{len(tends_ref) - 1} "
                f"tendencies: {', '.join(tends_fail)}."
            )

        diags_fail = validate(diags, diags_ref)
        if len(diags_fail) == 0:
            print("All diagnostics have been successfully validated. HOORAY!")
        else:
            print(
                f"Validation failed for {len(diags_fail)}/{len(diags_ref) - 1} "
                f"diagnostics: {', '.join(diags_fail)}."
            )


@click.command()
@click.option("--backend", type=str, default=None)
@click.option("--enable-checks/--disable-checks", is_flag=True, type=bool, default=False)
@click.option("--enable-validation/--disable-validation", is_flag=True, type=bool, default=True)
@click.option("--num-runs", type=int, default=None)
@click.option("--nx", type=int, default=None)
def main(
    backend: Optional[str],
    enable_checks: bool,
    enable_validation: bool,
    num_runs: Optional[int],
    nx: Optional[int],
) -> None:
    config = (
        default_config.with_backend(backend)
        .with_checks(enable_checks)
        .with_validation(enable_validation)
        .with_num_runs(num_runs)
        .with_nx(nx)
    )
    core(config)


if __name__ == "__main__":
    main()
