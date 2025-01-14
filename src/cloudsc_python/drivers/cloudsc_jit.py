# -*- coding: utf-8 -*-

# (C) Copyright 2018- ECMWF.

# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

import click
from pathlib import Path
import os
import sys
from ctypes import CDLL

from loki import config as loki_config
from loki.build import Builder, jit_compile_lib, Obj, Lib


@click.command()
@click.option(
    '--ngptot', default=100,
    help='Total number of columns to use for benchamrking'
)
@click.option(
    '--nproma', default=32,
    help='Number of columns per block (NPROMA)'
)
def main(ngptot, nproma):
    """
    Run a JIT-compiled version of CLOUDSC and validate against reference data
    """
    loki_config['log-level'] = 'debug'

    here = Path(__file__).parent.absolute()
    cloudsc_root = here.parent.parent.parent
    cloudsc_ftn = here.parent/'src/fortran'
    input_path = cloudsc_root/'config-files/input.h5'
    reference_path = cloudsc_root/'config-files/reference.h5'

    # The following assumes for now that Path.cwd() is dwarf-cloudsc/build
    build_dir = Path.cwd()
    tmp_path = build_dir/'jit'
    tmp_path.mkdir(exist_ok=True)

    builder = Builder(source_dirs=cloudsc_ftn, build_dir=tmp_path, workers=1)
    objs = [
        Obj(source_path=cloudsc_ftn/'cloudsc.F90'),
        Obj(source_path=cloudsc_ftn/'yomcst.F90'),
        Obj(source_path=cloudsc_ftn/'yoethf.F90'),
        Obj(source_path=cloudsc_ftn/'yoecldp.F90'),
    ]
    include_dirs = [
        tmp_path,
        build_dir/'cloudsc-dwarf/src/common/module',
        cloudsc_root/'src/common/module',
        cloudsc_root/'src/common/include',
    ]
    lib = Lib(name='cloudsc_f', objs=objs, shared=False)
    lib.build(builder=builder, include_dirs=include_dirs)
    sources = [
        cloudsc_ftn/'cloudsc.F90',
        cloudsc_ftn/'yomcst.F90',
        cloudsc_ftn/'yoethf.F90',
        cloudsc_ftn/'yoecldp.F90',
    ]

    # Pre-load the dependency library
    CDLL(build_dir/'lib64/libcloudsc-common-lib.so')
        
    cloudsc_fc = lib.wrap(
        modname='cloudsc_fc', sources=sources, builder=builder,
        libs=['cloudsc-common-lib'], lib_dirs=[build_dir/'lib64/'],
        kind_map=cloudsc_ftn/'kind_map'
    )
    cloudsc = cloudsc_fc.cloudsc_mod.cloudsc

    from IPython import embed; embed()
