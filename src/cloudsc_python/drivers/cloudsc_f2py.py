# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import click
from pathlib import Path


def loki_generate_kernel(source_path, out_path, include_dir=None):
    from loki import Sourcefile, flatten
    from loki.transform import FortranPythonTransformation

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
        includes=source_dir/'include', preprocess=True
    )
    f2py.apply(kernel, role='kernel', path=out_path)


def run_cloudsc_kernel(nthreads, ngptot, nproma, input_path, reference_path):
    from cloudscf2py import (
        load_input_fields, load_input_parameters, load_reference_fields,
        cloudsc_py
    )

    fields = load_input_fields(path=input_path)

    yrecldp, yrmcst, yrethf, yrephli, yrecld = load_input_parameters(path=input_path)

    cloudsc_args = {

    }
    cloudsc_args.update( (k.lower(), v) for k, v in fields.items() )

    cloudsc_py(
        kidia=1, kfdia=ngptot, **cloudsc_args,
        yrecldp=yrecldp, ydcst=yrmcst, ydthf=yrethf,
    )

    reference = load_reference_fields(path=reference_path)


@click.command()
@click.option(
    '--nthreads', default=1,
    help='Number of OpenMP threads to use'
)
@click.option(
    '--ngptot', default=100,
    help='Total number of columns to use for benchamrking'
)
@click.option(
    '--nproma', default=32,
    help='Number of columns per block (NPROMA)'
)
@click.option(
    '--generate/--no-generate', default=False,
    help='(Re)generate kernel via Loki-Fortran-Python transform'
)
def dwarf_cloudsc(nthreads, ngptot, nproma, generate):
    """
    Run a Python version of CLOUDSC and validate against reference data
    """

    here = Path(__file__).parent.absolute()
    cloudsc_root = here.parent.parent.parent
    cloudsc_f2py = here.parent/'src/cloudscf2py'
    input_path = cloudsc_root/'config-files/input.h5'
    reference_path = cloudsc_root/'config-files/reference.h5'

    if generate:
        loki_generate_kernel(
            source_path=cloudsc_f2py/'cloudsc.F90', out_path=cloudsc_f2py,
            include_dir=cloudsc_root/'src/common/include'
        )

    run_cloudsc_kernel(
        nthreads, ngptot, nproma, input_path=input_path, reference_path=reference_path
    )
        

if __name__ == "__main__":
    dwarf_cloudsc()
