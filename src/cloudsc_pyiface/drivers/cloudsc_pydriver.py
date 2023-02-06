"""
Driver that:
- loads input variables and parameters from .h5 file,
- invokes Fortran kernel computation,
- validates against reference results read from another .h5 file.
"""
import sys
import os
from pathlib import Path
from operator import itemgetter
from cloudsc_data import define_fortran_fields
from cloudsc_data import load_input_fortran_fields, load_input_parameters
from cloudsc_data import load_reference_fields
from cloudsc_data import convert_fortran_output_to_python
from cloudsc_data import cloudsc_validate
from importlib import import_module
here = os.getcwd()
cldir = here + '/../../cloudsc-dwarf/src/cloudsc_pyiface'
if cldir not in sys.path:
    sys.path.append(cldir)
clsc = import_module('cloudsc')

NPROMA=100
NUMOMP=1
NLEV=137
NGPTOT=100
NGPTOTG=100
NBLOCKS=1
PTSPHY=3600.

clsfields=define_fortran_fields(NPROMA,NLEV,NBLOCKS)

for fieldname in clsfields.keys():
    locals()[fieldname]=itemgetter(fieldname)(clsfields)



rootpath = Path(__file__).resolve().parents[3]
input_path = rootpath/'config-files/input.h5'


# Get referennce solution fields from file
ref_path = rootpath/'config-files/reference.h5'
ref_fields = load_reference_fields(path=ref_path)


NCLV = 5      # number of microphysics variables

load_input_parameters(input_path, ydecldp, ydephli, ydomcst, ydoethf)

input_fort_fields = load_input_fortran_fields(input_path,NPROMA,NLEV,NBLOCKS,clsfields)

for fieldname in input_fort_fields.keys():
    locals()[fieldname]=input_fort_fields[fieldname]

clsc.cloudsc_driver_mod.cloudsc_driver(
                         NUMOMP, NPROMA, NLEV, NGPTOT, NGPTOTG, NCLV,
                         kfldx, PTSPHY,
                         pt, pq,
                         buffer_tmp, buffer_loc,
                         pvfa, pvfl, pvfi, pdyna, pdynl, pdyni,
                         phrsw, phrlw,
                         pvervel, pap, paph,
                         plsm, ldcum, ktype,
                         plu, plude, psnde, pmfu, pmfd,
                         pa, pclv, psupsat,
                         plcrit_aer, picrit_aer, pre_ice,
                         pccn, pnice,
                         pcovptot, prainfrac_toprfz,
                         pfsqlf,   pfsqif ,  pfcqnng,  pfcqlng,
                         pfsqrf,   pfsqsf ,  pfcqrng,  pfcqsng,
                         pfsqltur, pfsqitur,
                         pfplsl, pfplsn, pfhpsl, pfhpsn,
                         ydomcst, ydoethf, ydecldp)

output_fields = convert_fortran_output_to_python (NPROMA,NLEV,NBLOCKS,input_fort_fields)

print ("Python-side validation:")
cloudsc_validate(output_fields, ref_fields)