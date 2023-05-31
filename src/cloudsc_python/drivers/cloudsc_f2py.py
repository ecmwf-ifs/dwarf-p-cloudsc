# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import click
from pathlib import Path
import numpy as np


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


def cloudsc_validate(fields, ref_fields, kidia, kfdia):
    _field_names = [
        'plude', 'pcovptot', 'prainfrac_toprfz', 'pfsqlf', 'pfsqif',
        'pfcqlng', 'pfcqnng', 'pfsqrf', 'pfsqsf', 'pfcqrng', 'pfcqsng',
        'pfsqltur', 'pfsqitur', 'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn',
        'tendency_loc_a', 'tendency_loc_q', 'tendency_loc_t',
        'tendency_loc_cld'
    ]
    ngptot = kfdia - kidia + 1

    print(
        "             Variable Dim             MinValue             MaxValue"
        "            AbsMaxErr         AvgAbsErr/GP          MaxRelErr-%"
    )
    for name in _field_names:
        if len(fields[name].shape) == 1:
            f = fields[name][kidia-1:kfdia]
            ref = ref_fields[name][kidia-1:kfdia]
        elif len(fields[name].shape) == 2:
            f = fields[name][:,kidia-1:kfdia]
            ref = ref_fields[name][:,kidia-1:kfdia]
        elif len(fields[name].shape) == 3:
            f = fields[name][:,:,kidia-1:kfdia]
            ref = ref_fields[name][:,:,kidia-1:kfdia]
        else:
            f = fields[name]
            ref = ref_fields[name]
        zsum = np.sum(np.absolute(ref))
        zerrsum = np.sum(np.absolute(f - ref))
        zeps = np.finfo(np.float64).eps
        print(
            ' {fname:>20}     {fmin:20.13e}  {fmax:20.13e}  {absmax:20.13e} '
            ' {absavg:20.13e}  {maxrel:20.13e}'.format(
                fname=name.upper(), fmin=f.min(), fmax=f.max(),
                absmax=np.absolute(f - ref).max(),
                absavg=np.sum(np.absolute(f - ref)) / ngptot,
                maxrel=0.0 if zerrsum < zeps else (zerrsum/(1.0+zsum) if zsum < zeps else zerrsum/zsum)
            )
        )


def run_cloudsc_kernel(ngptot, nproma, input_path, reference_path):
    from cloudscf2py import (
        load_input_fields, load_input_parameters, load_reference_fields,
        cloudsc_py
    )

    fields = load_input_fields(path=input_path, ngptot=ngptot)

    yrecldp, yrmcst, yrethf, yrephli, yrecld = load_input_parameters(path=input_path)

    cloudsc_args = {k.lower(): v for k, v in fields.items()}

    # We process only one block for now, all in one go
    cloudsc_args['klon'] = ngptot

    cloudsc_py(
        kidia=1, kfdia=ngptot, **cloudsc_args,
        yrecldp=yrecldp, ydcst=yrmcst, ydthf=yrethf,
    )

    # Validate the output fields against reference data
    reference = load_reference_fields(path=reference_path, ngptot=ngptot)
    cloudsc_validate(cloudsc_args, reference, kidia=1, kfdia=ngptot)


@click.command()
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
def main(ngptot, nproma, generate):
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
        ngptot, nproma, input_path=input_path, reference_path=reference_path
    )
        

if __name__ == "__main__":
    dwarf_cloudsc()
