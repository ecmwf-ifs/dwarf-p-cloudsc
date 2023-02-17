#!/usr/bin/env python3.11
import sys,os
from pathlib import Path
from operator import itemgetter
from cloudscpytools import cloudsc_data
from importlib import import_module
import click

@click.command()
@click.option(
    "--numomp",
    type=int,
    default=1,
    help="Python driver to execute IFS physics kernel (CLOUDSC)."
#    "\n\nOptions: numomp, ngptot, nproma",
#    "\n\nNumber of OPenMP threads used for the benchmark. Default: 1",
)
@click.option(
    "--ngptot",
    type=int,
    default=100,
    help="Total number of grid points (NGPTOT) used for the benchmark. Default: 100",
)
@click.option(
    "--nproma",
    type=int,
    default=100,
    help="Block sizes (NPROMA) used for the benchmark. Default: 100",
)
def main(
    numomp: int,
    ngptot: int,
    nproma: int
) -> None:
    """
    Driver that:
    - loads input variables and parameters from .h5 file,
    - invokes Fortran kernel computation,
    - validates against reference results read from another .h5 file.
    """
    
    here = os.getcwd()
    cldir = here + '/../../cloudsc-dwarf/src/cloudsc_pyiface'
    if cldir not in sys.path:
        sys.path.append(cldir)
    clsc = import_module('cloudsc')
    
    NPROMA=nproma
    NUMOMP=numomp
    NLEV=137
    NGPTOT=ngptot
    NGPTOTG=ngptot
    NBLOCKS= int( (ngptot / nproma) + min(ngptot % nproma, 1) ) 
    PTSPHY=3600.
    
    clsfields=cloudsc_data.define_fortran_fields(NPROMA,NLEV,NBLOCKS)
    
    
    rootpath = Path(__file__).resolve().parents[3]
    input_path = rootpath/'config-files/input.h5'
    
    
    # Get referennce solution fields from file
    ref_path = rootpath/'config-files/reference.h5'
    ref_fields = cloudsc_data.load_reference_fields(path=ref_path)
    
    
    NCLV = 5      # number of microphysics variables
    cloudsc_data.load_input_parameters(input_path, itemgetter('ydecldp')(clsfields),
                                                   itemgetter('ydephli')(clsfields), 
                                                   itemgetter('ydomcst')(clsfields),
                                                   itemgetter('ydoethf')(clsfields))
    
    input_fort_fields = cloudsc_data.load_input_fortran_fields(input_path,NPROMA,NLEV,NBLOCKS,clsfields)
    
    clsc.cloudsc_driver_mod.cloudsc_driver(
                             NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NCLV,
                             itemgetter('kfldx')(input_fort_fields), 
                             itemgetter('PTSPHY')(input_fort_fields),
                             itemgetter('pt')(input_fort_fields),
                             itemgetter('pq')(input_fort_fields),
                             itemgetter('buffer_tmp')(clsfields),
                             itemgetter('buffer_loc')(clsfields),
                             itemgetter('pvfa')(input_fort_fields), 
                             itemgetter('pvfl')(input_fort_fields), 
                             itemgetter('pvfi')(input_fort_fields), 
                             itemgetter('pdyna')(input_fort_fields), 
                             itemgetter('pdynl')(input_fort_fields), 
                             itemgetter('pdyni')(input_fort_fields),
                             itemgetter('phrsw')(input_fort_fields), 
                             itemgetter('phrlw')(input_fort_fields),
                             itemgetter('pvervel')(input_fort_fields), 
                             itemgetter('pap')(input_fort_fields), 
                             itemgetter('paph')(input_fort_fields),
                             itemgetter('plsm')(input_fort_fields), 
                             itemgetter('ldcum')(input_fort_fields), 
                             itemgetter('ktype')(input_fort_fields),
                             itemgetter('plu')(input_fort_fields), 
                             itemgetter('plude')(input_fort_fields), 
                             itemgetter('psnde')(input_fort_fields), 
                             itemgetter('pmfu')(input_fort_fields), 
                             itemgetter('pmfd')(input_fort_fields),
                             itemgetter('pa')(input_fort_fields), 
                             itemgetter('pclv')(input_fort_fields), 
                             itemgetter('psupsat')(input_fort_fields),
                             itemgetter('plcrit_aer')(input_fort_fields), 
                             itemgetter('picrit_aer')(input_fort_fields), 
                             itemgetter('pre_ice')(input_fort_fields),
                             itemgetter('pccn')(input_fort_fields), 
                             itemgetter('pnice')(input_fort_fields),
                             itemgetter('pcovptot')(clsfields), 
                             itemgetter('prainfrac_toprfz')(clsfields),
                             itemgetter('pfsqlf')(clsfields),   
                             itemgetter('pfsqif')(clsfields),  
                             itemgetter('pfcqnng')(clsfields),  
                             itemgetter('pfcqlng')(clsfields),
                             itemgetter('pfsqrf')(clsfields),   
                             itemgetter('pfsqsf')(clsfields) ,  
                             itemgetter('pfcqrng')(clsfields),  
                             itemgetter('pfcqsng')(clsfields),
                             itemgetter('pfsqltur')(clsfields), 
                             itemgetter('pfsqitur')(clsfields),
                             itemgetter('pfplsl')(clsfields), 
                             itemgetter('pfplsn')(clsfields), 
                             itemgetter('pfhpsl')(clsfields), 
                             itemgetter('pfhpsn')(clsfields),
                             itemgetter('ydomcst')(clsfields), 
                             itemgetter('ydoethf')(clsfields), 
                             itemgetter('ydecldp')(clsfields))
    
    output_fields = cloudsc_data.convert_fortran_output_to_python (NPROMA,NLEV,NBLOCKS,clsfields)
    
    print ("Python-side validation:")
    cloudsc_data.cloudsc_validate(output_fields, ref_fields)


if __name__ == "__main__":
    main()
