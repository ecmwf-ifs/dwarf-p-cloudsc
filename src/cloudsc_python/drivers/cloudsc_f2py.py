# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import click
from pathlib import Path


def loki_generate_kernel(source_path, out_path, include_dir=None, log_level='perf'):
    from loki import Sourcefile, flatten, config as loki_config
    from loki.transformations import FortranPythonTransformation

    loki_config['log-level'] = log_level

    source_dir = source_path.parent
    headers = ['yoethf.F90', 'yoecldp.F90', 'yomcst.F90']
    definitions = flatten(
        Sourcefile.from_file(source_dir/header).modules for header in headers
    )

    f2py = FortranPythonTransformation(
        with_dace=False, suffix='_py', invert_indices=True
    )

    # Parse original driver and kernel routine, and enrich the driver
    kernel = Sourcefile.from_file(
        source_path, definitions=definitions,
        includes=include_dir, preprocess=True
    )
    f2py.apply(kernel['cloudsc'], role='kernel', path=out_path)


def run_cloudsc_kernel(ngptot, nproma, input_path, reference_path):
    from cloudscf2py import (
        load_input_fields, load_input_parameters, load_reference_fields,
        validate, cloudsc_py
    )

    fields = load_input_fields(path=input_path, ngptot=ngptot)

    yrecldp, yrmcst, yrethf = load_input_parameters(path=input_path)

    cloudsc_args = {k.lower(): v for k, v in fields.items()}

    # We process only one block for now, all in one go
    cloudsc_args['klon'] = ngptot

    cloudsc_py(
        kidia=1, kfdia=ngptot, **cloudsc_args,
        yrecldp=yrecldp, ydcst=yrmcst, ydthf=yrethf,
    )

    # Validate the output fields against reference data
    reference = load_reference_fields(path=reference_path, ngptot=ngptot)
    validate(cloudsc_args, reference, kidia=1, kfdia=ngptot)


@click.command()
@click.option(
    '--ngptot', default=100,
    help='Total number of columns to use for benchmarking'
)
@click.option(
    '--nproma', default=32,
    help='Number of columns per block (NPROMA)'
)
@click.option(
    '--generate/--no-generate', default=False,
    help='(Re)generate kernel via Loki-Fortran-Python transform'
)
@click.option(
    '--log-level', '-l', default='perf',
    help='Log-level to set for Loki when regenerating kernel'
)
def main(ngptot, nproma, generate, log_level):
    """
    Run a Python version of CLOUDSC and validate against reference data
    """

    here = Path(__file__).parent.absolute()
    cloudsc_root = here.parent.parent.parent
    cloudsc_python = here.parent
    cloudsc_f2py = here.parent/'src/cloudscf2py'
    input_path = cloudsc_root/'config-files/input.h5'
    reference_path = cloudsc_root/'config-files/reference.h5'

    if generate:
        loki_generate_kernel(
            source_path=cloudsc_python/'src/fortran/cloudsc.F90',
            include_dir=cloudsc_python/'src/include',
            out_path=cloudsc_f2py, log_level=log_level
        )

    run_cloudsc_kernel(
        ngptot, nproma, input_path=input_path, reference_path=reference_path
    )
        

if __name__ == "__main__":
    dwarf_cloudsc()
