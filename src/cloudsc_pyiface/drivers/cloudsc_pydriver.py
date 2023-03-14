#!/usr/bin/env python3
"""
Driver that executes Fortran implementation of the CLOUDSC dwarf using f90wrap/f2py
"""
import sys
import os
from pathlib import Path
from importlib import import_module
import click

from cloudscpytools import cloudsc_data

@click.command()
@click.option(
    "--numomp", type=int, default=1,
    help="Number of OpenMP threads used for the benchmark. Default: 1",
)
@click.option(
    "--ngptot", type=int, default=100,
    help="Total number of grid points (NGPTOT) used for the benchmark. Default: 100",
)
@click.option(
    "--nproma", type=int, default=100,
    help="Block sizes (NPROMA) used for the benchmark. Default: 100",
)
def main(numomp: int, ngptot: int, nproma: int) -> None:
    """
    Python driver to execute IFS physics kernel (CLOUDSC).

    Performs the following tasks:
    - loads input variables and parameters from .h5 file,
    - invokes Fortran kernel computation,
    - validates against reference results read from another .h5 file.
    """
    here = os.getcwd()

    # Dynamically loading Python-wrapped CLOUDSC Fortran module
    cldir = here + '/../../cloudsc-dwarf/src/cloudsc_pyiface'
    if cldir not in sys.path:
        sys.path.append(cldir)
    clsc = import_module('cloudsc')

    # Defining common parameters
    nlev = 137
    ndim = 5
    ngptotg = ngptot
    nblocks = int( (ngptot / nproma) + min(ngptot % nproma, 1) )
    nclv = 5      # number of microphysics variables
    npars = dict(
        nlev=nlev, ngptot=ngptot, ngptotg=ngptotg, nproma=nproma,
        nblocks=nblocks, ndim=ndim, nclv=nclv
    )

    # Allocate temporary and output fields
    clsfields = cloudsc_data.define_fortran_fields(
        nproma=nproma, nlev=nlev, nblocks=nblocks
    )

    # Get reference solution fields from file
    rootpath = Path(__file__).resolve().parents[3]
    ref_path = rootpath/'config-files/reference.h5'
    ref_fields = cloudsc_data.load_reference_fields(path=ref_path, **npars)

    # Get input data fields from file
    input_path = rootpath/'config-files/input.h5'
    cloudsc_data.load_input_parameters(
        input_path, clsfields['ydecldp'], clsfields['ydephli'],
        clsfields['ydomcst'], clsfields['ydoethf']
    )
    input_fort_fields = cloudsc_data.load_input_fortran_fields(
        path=input_path, fields=clsfields, **npars
    )

    # Execute kernel via Python-wrapped, compiled Fortran driver
    clsc.cloudsc_driver_mod.cloudsc_driver(
        numomp, nproma, nlev, ngptot, ngptotg, nclv,
        input_fort_fields['kfldx'],
        input_fort_fields['PTSPHY'],
        input_fort_fields['pt'],
        input_fort_fields['pq'],
        clsfields['buffer_tmp'],
        clsfields['buffer_loc'],
        input_fort_fields['pvfa'],
        input_fort_fields['pvfl'],
        input_fort_fields['pvfi'],
        input_fort_fields['pdyna'],
        input_fort_fields['pdynl'],
        input_fort_fields['pdyni'],
        input_fort_fields['phrsw'],
        input_fort_fields['phrlw'],
        input_fort_fields['pvervel'],
        input_fort_fields['pap'],
        input_fort_fields['paph'],
        input_fort_fields['plsm'],
        input_fort_fields['ldcum'],
        input_fort_fields['ktype'],
        input_fort_fields['plu'],
        input_fort_fields['plude'],
        input_fort_fields['psnde'],
        input_fort_fields['pmfu'],
        input_fort_fields['pmfd'],
        input_fort_fields['pa'],
        input_fort_fields['pclv'],
        input_fort_fields['psupsat'],
        input_fort_fields['plcrit_aer'],
        input_fort_fields['picrit_aer'],
        input_fort_fields['pre_ice'],
        input_fort_fields['pccn'],
        input_fort_fields['pnice'],
        clsfields['pcovptot'],
        clsfields['prainfrac_toprfz'],
        clsfields['pfsqlf'],
        clsfields['pfsqif'],
        clsfields['pfcqnng'],
        clsfields['pfcqlng'],
        clsfields['pfsqrf'],
        clsfields['pfsqsf'],
        clsfields['pfcqrng'],
        clsfields['pfcqsng'],
        clsfields['pfsqltur'],
        clsfields['pfsqitur'],
        clsfields['pfplsl'],
        clsfields['pfplsn'],
        clsfields['pfhpsl'],
        clsfields['pfhpsn'],
        clsfields['ydomcst'],
        clsfields['ydoethf'],
        clsfields['ydecldp'],
    )
    output_fields = cloudsc_data.convert_fortran_output_to_python (clsfields, **npars)
    print ("Python-side validation:")
    cloudsc_data.cloudsc_validate(output_fields, ref_fields)


if __name__ == "__main__":
    main()
