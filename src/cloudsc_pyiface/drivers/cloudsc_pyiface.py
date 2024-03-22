#!/usr/bin/env python3
"""
Driver that executes Fortran implementation of the CLOUDSC dwarf using f90wrap/f2py
"""
from pathlib import Path
import click

from pyiface import cloudsc_data
from pyiface.dynload import load_module


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
@click.option(
    "--cloudsc-path", type=click.Path(exists=True), default=Path.cwd(),
    help="Path to the Python-wrapped and compiled CLOUDSC module",
)
@click.option(
    "--input-path", type=click.Path(exists=True), default=Path.cwd(),
    help="Path to input and reference files; by default './'",
)
def main(numomp: int, ngptot: int, nproma: int, cloudsc_path, input_path) -> None:
    """
    Python driver to execute IFS physics kernel (CLOUDSC).

    Performs the following tasks:
    - loads input variables and parameters from .h5 file,
    - invokes Fortran kernel computation,
    - validates against reference results read from another .h5 file.
    """

    cloudsc_path = Path(cloudsc_path)
    input_path = Path(input_path)

    # Dynamically load the Python-wrapped Fortran CLOUDSC module
    clsc = load_module(module='cloudsc', modpath=Path(cloudsc_path))

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
        nproma=nproma, nlev=nlev, nblocks=nblocks, clsc=clsc
    )

    # Get reference solution fields from file
    ref_fields = cloudsc_data.load_reference_fields(
        path=input_path/'reference.h5', clsc=clsc, **npars
    )

    # Get input data fields from file
    cloudsc_data.load_input_parameters(
        input_path/'input.h5', clsfields['ydecldp'], clsfields['ydephli'],
        clsfields['ydomcst'], clsfields['ydoethf']
    )
    input_fort_fields = cloudsc_data.load_input_fortran_fields(
        path=input_path/'input.h5', fields=clsfields, clsc=clsc, **npars
    )

    # Execute kernel via Python-wrapped, compiled Fortran driver
    clsc.cloudsc_driver_mod.cloudsc_driver(
        numomp, nproma, nlev, ngptot, ngptotg,
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
