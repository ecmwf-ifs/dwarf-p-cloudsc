"""
cloudsc_data module consist of utilities that:
- load variables serving as an input to a Fortran computational kernel;
- load physical parameters needed by a Fortran kernel;
- load reference results that will be compared with an output of Fortran computation;
- validates reference vs. computed fields;
- other, purely technical utilities.
"""
from collections import OrderedDict
import h5py
import numpy as np

NCLV = 5      # number of microphysics variables


def define_fortran_fields(nproma, nlev, nblocks, clsc):
    """
    define_fortran_fields returns:
    - zero NumPy arrays that will further be used as an output of Fortran kernel computation.
    - empty Fortran paramter datatypes that are created used constructors supplied by f90wrap.
    """

    fields = OrderedDict()

    argnames_nlev = [
        'pcovptot'
    ]

    argnames_nlevp = [
        'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn',
        'pfsqlf',   'pfsqif' ,  'pfcqnng',  'pfcqlng',
        'pfsqrf',   'pfsqsf' ,  'pfcqrng',  'pfcqsng',
        'pfsqltur', 'pfsqitur'
    ]

    argnames_buffer = [
        'buffer_loc','buffer_tmp'
    ]

    argnames_tend = [
        'tendency_loc_a','tendency_loc_t','tendency_loc_q',
    ]

    argnames_tend_cld = [
        'tendency_loc_cld'
    ]

    argnames_nproma = [
        'prainfrac_toprfz'
    ]

    for argname in argnames_nlev:
        fields[argname] = np.zeros(shape=(nproma,nlev  ,nblocks), order='F')

    for argname in argnames_nlevp:
        fields[argname] = np.zeros(shape=(nproma,nlev+1,nblocks), order='F')

    for argname in argnames_buffer:
        fields[argname] = np.zeros(shape=(nproma,nlev,3+NCLV,nblocks), order='F')

    for argname in argnames_tend:
        fields[argname] = np.zeros(shape=(nproma,nlev,nblocks), order='F')

    for argname in argnames_tend_cld:
        fields[argname] = np.zeros(shape=(nproma,nlev,NCLV,nblocks), order='F')


    for argname in argnames_nproma:
        fields[argname] = np.zeros(shape=(nproma,nblocks), order='F')

    fields['ydomcst']=clsc.yomcst.TOMCST()
    fields['ydoethf']=clsc.yoethf.TOETHF()
    fields['ydecldp']=clsc.yoecldp.TECLDP()
    fields['ydephli']=clsc.yoephli.TEPHLI()

    return fields


def field_c_to_fortran(dims, cfield, clsc=None, **kwargs):
    """
    field_c_to_fortran:
    1) transposes C array input to Fortran array
    2) rewrites Fortran linear array into block structure
    """

    # Transpose the C array (row-major)  into Fortran (column-major) data layout
    ffieldtmp = np.asfortranarray(np.transpose(np.ascontiguousarray(cfield)))

    return field_linear_to_block(dims, ffieldtmp, clsc=clsc, **kwargs)


def field_linear_to_block(dims, lfield, clsc=None, **kwargs):
    """
    Rewrites Fortran linear array into block structure
    """

    if not clsc:
        raise RuntimeError('[PyIface] Cannot expand field without CLOUDSC Fortran backend')

    # Pick array dimension arguments from keyword args
    nproma = kwargs.get('nproma', 32)
    nblocks = kwargs.get('nblocks', 1)
    ngptot = kwargs.get('ngptot', 100)
    nlon = kwargs.get('nlon', 100)
    ndim = kwargs.get('ndim', 1)
    nlev = dims[-2] #nparms['NLEV']
    ldims = len(dims)

    if lfield.dtype == "float64":
        if ldims == 2:
            b2field=np.asfortranarray(np.transpose(
                                      np.zeros(shape=dims, dtype="float64")))
            clsc.expand_mod.expand_r1(lfield, b2field, nlon, nproma, ngptot, nblocks)
            bfield=b2field
        elif ldims == 3:
            b3field=np.asfortranarray(np.transpose(np.zeros(shape=dims, dtype="float64")))
            clsc.expand_mod.expand_r2(lfield, b3field, nlon, nproma, nlev, ngptot, nblocks)
            bfield=b3field
        elif ldims == 4:
            b4field=np.asfortranarray(np.transpose(np.zeros(shape=dims, dtype="float64")))
            clsc.expand_mod.expand_r3(lfield, b4field, nlon=nlon, nproma=nproma, nlev=nlev,
                                      ndim=ndim, ngptot=ngptot, nblocks=nblocks)
            bfield=b4field
        else:
            print ("Wrong float ldim")
    elif lfield.dtype == "bool":
        # Workaround - using type int32, otherwise complains about type disagreement at runtime
        bfield=np.asfortranarray(np.transpose(np.zeros(shape=dims, dtype='int32')))
        if ldims == 2:
            tlfield=lfield.astype('int32')
            clsc.expand_mod.expand_l1(tlfield, bfield,  nlon, nproma, ngptot, nblocks)
        else:
            print ("Wrong bool ldim")
    elif lfield.dtype == "int32":
        bfield=np.asfortranarray(np.transpose(np.zeros(shape=dims, dtype='int32')))
        if ldims == 2:
            clsc.expand_mod.expand_i1(lfield, bfield,  nlon, nproma, ngptot, nblocks)
        else:
            print ("Wrong int ldim")
    else:
        print ("Wrong dtype")
    return bfield

def load_input_fortran_fields(path, fields, clsc=None, **kwargs):
    """
    load_input_fortran_fields returns:
    - set of variables needed to initiate computation of the Fortran kernel.
    """

    if not clsc:
        raise RuntimeError('[PyIface] Cannot load input fields without CLOUDSC Fortran backend')
    
    nproma = kwargs['nproma']
    nlev = kwargs['nlev']
    nblocks = kwargs['nblocks']
    argnames_nlev = [
        'pt', 'pq',
        'pvfa', 'pvfl', 'pvfi', 'pdyna', 'pdynl', 'pdyni',
        'phrsw', 'phrlw','pvervel','pap','plu','plude',
        'psnde', 'pmfu', 'pmfd',
        'pa', 'psupsat',
        'plcrit_aer','picrit_aer','pre_ice',
        'pccn', 'pnice'
    ]
    argnames_nlevp = [
        'paph'
    ]

    argnames_withnclv= [
        'pclv','tendency_tmp_cld'
    ]

    argnames_tend = [
        'tendency_tmp_t','tendency_tmp_q','tendency_tmp_a'
    ]

    argnames_scalar = [
         'kfldx'
    ]

    argnames_nproma = [
        'plsm', 'ldcum', 'ktype'
    ]


    with h5py.File(path, 'r') as f:
        fields['KLON'] = f['KLON'][0]
        fields['KLEV'] = f['KLEV'][0]
        fields['PTSPHY'] = f['PTSPHY'][0]
        kwargs['nlon'] = fields['KLON']

        for argname in argnames_nlev:
            fields[argname] = field_c_to_fortran((nblocks,nlev,nproma),
                                                 f[argname.upper()], clsc=clsc, **kwargs)

        for argname in argnames_nlevp:
            fields[argname] = field_c_to_fortran((nblocks,nlev+1,nproma),
                                                 f[argname.upper()], clsc=clsc, **kwargs)

        for argname in argnames_withnclv:
            fields[argname] = field_c_to_fortran((nblocks,NCLV,nlev,nproma),
                                                 f[argname.upper()], clsc=clsc, **kwargs)

        for argname in argnames_tend:
            fields[argname] = field_c_to_fortran((nblocks,nlev,nproma),
                                                 f[argname.upper()], clsc=clsc, **kwargs)

        for argname in argnames_nproma:
            fields[argname] = field_c_to_fortran((nblocks,nproma),
                                                 f[argname.upper()], clsc=clsc, **kwargs)

        for argname in argnames_scalar:
            fields[argname] = f[argname.upper()][0]

    pack_buffer_using_tendencies(fields['buffer_tmp'      ],
                                 fields['tendency_tmp_a'  ],
                                 fields['tendency_tmp_t'  ],
                                 fields['tendency_tmp_q'  ],
                                 fields['tendency_tmp_cld'])
    return fields

def pack_buffer_using_tendencies(buffervar,tendency_a,tendency_t,tendency_q,tendency_cld):
    """
    pack_buffer_using_tendencies serves as a packager of a single-variable
    (that may consist of multiple fields, e.g. moist species)
    tendencies into a continous buffer
    """
    buffervar[:,:,0       ,:]=tendency_t  [:,:,:]
    buffervar[:,:,1       ,:]=tendency_a  [:,:,:]
    buffervar[:,:,2       ,:]=tendency_q  [:,:,:]
    buffervar[:,:,3:3+NCLV-1,:]=tendency_cld[:,:,0:NCLV-1,:]

def unpack_buffer_to_tendencies(buffervar,tendency_a,tendency_t,tendency_q,tendency_cld):
    """
    unpack_buffer_to_tendencies continuous unpacks buffer into a set of a single-variable
    (that may consist of multiple fields, e.g. moist species) tendencies.
    """
    tendency_t  [:,:,:]=buffervar[:,:,0       ,:]
    tendency_a  [:,:,:]=buffervar[:,:,1       ,:]
    tendency_q  [:,:,:]=buffervar[:,:,2       ,:]
    tendency_cld[:,:,0:NCLV-1,:]=buffervar[:,:,3:3+NCLV-1,:]

def load_input_parameters(path,yrecldp,yrephli,yrmcst,yrethf):
    """
    load_input_parameters returns:
    - four parameter datatypes that are filled using names read from the reference .h5 file
    """
    with h5py.File(path, 'r') as f:
        tecldp_keys = [k for k in f.keys() if 'YRECLDP' in k]
        for k in tecldp_keys:
            attrkey = k.replace('YRECLDP_', '').lower()
            setattr(yrecldp, attrkey, f[k][0])
        yrecldp.ncldql = 1
        yrecldp.ncldqi = 2
        yrecldp.ncldqr = 3
        yrecldp.ncldqs = 4
        yrecldp.ncldqv = 5

        tephli_keys = [k for k in f.keys() if 'YREPHLI' in k]
        for k in tephli_keys:
            attrkey = k.replace('YREPHLI_', '').lower()
            setattr(yrephli, attrkey, f[k][0])

        tomcst_keys = ['RG', 'RD', 'RCPD', 'RETV', 'RLVTT', 'RLSTT', 'RLMLT', 'RTT', 'RV' ]
        for k in tomcst_keys:
            attrkey = k.lower()
            setattr(yrmcst, attrkey, f[k][0])

        toethf_keys = ['R2ES', 'R3LES', 'R3IES', 'R4LES', 'R4IES', 'R5LES', 'R5IES',
                       'R5ALVCP', 'R5ALSCP', 'RALVDCP', 'RALSDCP', 'RALFDCP',
                       'RTWAT', 'RTICE', 'RTICECU', 'RTWAT_RTICE_R', 'RTWAT_RTICECU_R',
                       'RKOOP1', 'RKOOP2'  ]

        for k in toethf_keys:
            attrkey = k.lower()
            setattr(yrethf, attrkey, f[k][0])

def convert_fortran_output_to_python (input_fields, **kwargs):
    """
    convert_fortran_output_to_python converts Fortran-format fields that are to be compared to
    reference results into a Python format.
    """
    nproma = kwargs['nproma']
    nlev = kwargs['nlev']
    nblocks = kwargs['nblocks']

    fields = OrderedDict()
    argnames_nlev = [
        'plude', 'pcovptot'
    ]

    argnames_nlevp = [
        'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn',
        'pfsqlf',   'pfsqif' ,  'pfcqnng',  'pfcqlng',
        'pfsqrf',   'pfsqsf' ,  'pfcqrng',  'pfcqsng',
        'pfsqltur', 'pfsqitur'
    ]

    argnames_nproma = [
        'prainfrac_toprfz'
    ]

    argnames_tend = [
        'tendency_loc_a','tendency_loc_t','tendency_loc_q'
    ]

    argnames_tend_cld = [
        'tendency_loc_cld'
    ]

    for argname in argnames_nlev:
        fields[argname] = input_fields[argname]

    for argname in argnames_nlevp:
        fields[argname] = input_fields[argname]

    for argname in argnames_nproma:
        fields[argname] = input_fields[argname]

    for argname in argnames_tend:
        fields[argname] = np.zeros(shape=(nproma,nlev,nblocks), order='F')

    for argname in argnames_tend_cld:
        fields[argname] = np.zeros(shape=(nproma,nlev,NCLV,nblocks), order='F')


    unpack_buffer_to_tendencies(input_fields ['buffer_loc'],
                                      fields ['tendency_loc_a'],
                                      fields ['tendency_loc_t'],
                                      fields ['tendency_loc_q'],
                                      fields ['tendency_loc_cld'])

    return fields

def load_reference_fields (path, clsc=None, **kwargs):
    """
    load_reference_fields loads reference results of Fortran computation from the .h5 file
    """

    if not clsc:
        raise RuntimeError('[PyIface] Cannot load reference fields without CLOUDSC Fortran backend')

    nproma = kwargs['nproma']
    nlev = kwargs['nlev']
    nblocks = kwargs['nblocks']
    fields = OrderedDict()

    argnames_nlev = [
        'plude', 'pcovptot'
    ]

    argnames_nproma = [
         'prainfrac_toprfz'
    ]

    argnames_nlevp = [
        'pfsqlf',   'pfsqif' ,  'pfcqnng',  'pfcqlng',
        'pfsqrf',   'pfsqsf' ,  'pfcqrng',  'pfcqsng',
        'pfsqltur', 'pfsqitur' ,
        'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn'
    ]

    argnames_tend = [
        'tendency_loc_a','tendency_loc_t','tendency_loc_q',
    ]

    argnames_tend_cld = [
        'tendency_loc_cld'
    ]

    with h5py.File(path, 'r') as f:
        fields['KLON'] = f['KLON'][0]
        fields['KLEV'] = f['KLEV'][0]
        kwargs['nlon'] = fields['KLON']

        for argname in argnames_nlev:
            fields[argname] = field_c_to_fortran((nblocks,nlev,nproma),
                                                 f[argname.upper()], clsc=clsc, **kwargs)

        for argname in argnames_nlevp:
            fields[argname] = field_c_to_fortran((nblocks,nlev+1,nproma),
                                                 f[argname.upper()], clsc=clsc, **kwargs)

        for argname in argnames_nproma:
            fields[argname] = field_c_to_fortran((nblocks,nproma),
                                                 f[argname.upper()], clsc=clsc, **kwargs)

        for argname in argnames_tend:
            fields[argname] = field_c_to_fortran((nblocks,nlev,nproma),
                                                 f[argname.upper()], clsc=clsc, **kwargs)

        for argname in argnames_tend_cld:
            fields[argname] = field_c_to_fortran((nblocks,NCLV,nlev,nproma),
                                                 f[argname.upper()], clsc=clsc, **kwargs)

    return fields

def cloudsc_validate(fields, ref_fields):
    """
    cloudsc_validate compares computed output of a Fortran kernel with reference results
    previously read from the .h5 file.
    """
    # List of refencece fields names in order
    _field_names = [
        'plude', 'pcovptot','prainfrac_toprfz',
        'pfsqlf',   'pfsqif' ,  'pfcqnng',  'pfcqlng',
        'pfsqrf',   'pfsqsf' ,  'pfcqrng',  'pfcqsng',
        'pfsqltur', 'pfsqitur' ,
        'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn',
        'tendency_loc_a', 'tendency_loc_q', 'tendency_loc_t', 'tendency_loc_cld'
    ]
    kidia = 1
    kfdia = 100
    ngptot = kfdia - kidia + 1

    print("             Variable Dim             MinValue             MaxValue\
            AbsMaxErr         AvgAbsErr/GP          MaxRelErr-%")
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
        print(' {fname:>20}     {fmin:20.13e}  {fmax:20.13e}  {absmax:20.13e} '\
              ' {absavg:20.13e}  {maxrel:20.13e}'.format(
                  fname=name.upper(), fmin=f.min(), fmax=f.max(),
                  absmax=np.absolute(f - ref).max(),
                  absavg=np.sum(np.absolute(f - ref)) / ngptot,
                  maxrel=0.0 if zerrsum < zeps else (zerrsum/(1.0+zsum)
                             if    zsum < zeps else zerrsum/zsum)
              )
        )
