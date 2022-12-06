import h5py
import numpy as np
import cloudsc2 as clsc
from pathlib import Path
from collections import OrderedDict


NCLV = 5      # number of microphysics variables


def define_fortran_fields(nproma,nlev,nblocks):
    """

    """
    fields = OrderedDict()

    argnames_nlev = [
        'pt', 'pq', 'pap', 'plu', 'plude', 'pmfu', 'pmfd',
        'pa', 'psupsat', 'pcovptot'
    ]

    argnames_nlevp = [
        'paph', 'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn'
    ]

    argnames_withnclv= [
        'pclv','tendency_loc_cld'
    ]

    argnames_buffer = [
        'buffer_cml', 'buffer_loc'
    ]
 
    argnames_tend = [
        'tendency_loc_a','tendency_loc_t','tendency_loc_q'
    ]

    for argname in argnames_nlev:
        fields[argname] = np.zeros(shape=(nproma,nlev  ,nblocks), order='F')

    for argname in argnames_nlevp:
        fields[argname] = np.zeros(shape=(nproma,nlev+1,nblocks), order='F')

    for argname in argnames_withnclv:
        fields[argname] = np.zeros(shape=(nproma,nlev,NCLV,nblocks), order='F')

    for argname in argnames_buffer:
        fields[argname] = np.zeros(shape=(nproma,nlev,3+NCLV,nblocks), order='F')

    for argname in argnames_tend:
        fields[argname] = np.zeros(shape=(nproma,nlev,nblocks), order='F')

    fields['ydomcst']=clsc.yomcst.TOMCST()
    fields['ydoethf']=clsc.yoethf.TOETHF()
    fields['ydecld' ]=clsc.yoecld.TECLD()
    fields['ydecldp']=clsc.yoecldp.TECLDP()
    fields['ydephli']=clsc.yoephli.TEPHLI()
    fields['ydphnc' ]=clsc.yophnc.TPHNC()

    return fields

def load_input_fortran_fields(path, nproma, nlev, nblocks,  transpose=False):
    """

    """
    fields = OrderedDict()

    argnames_nlev = [
        'pt', 'pq', 'pap', 'plu', 'plude', 'pmfu', 'pmfd',
        'pa', 'psupsat', 'pcovptot'
    ]

    argnames_nlevp = [
        'paph'
    ]

    argnames_withnclv= [
        'pclv','tendency_cml_cld'
    ]

    argnames_tend = [
        'tendency_cml_t','tendency_cml_q'
    ]

    argnames = [
        'pt', 'pq', 'pap', 'paph', 'plu', 'plude', 'pmfu', 'pmfd',
        'pa', 'pclv', 'psupsat', 'tendency_cml_t', 'tendency_cml_q',
        'tendency_cml_cld'
    ]

    with h5py.File(path, 'r') as f:
        fields['KLON'] = f['KLON'][0]
        fields['KLEV'] = f['KLEV'][0]
        fields['PTSPHY'] = f['PTSPHY'][0]

        klon = fields['KLON']
        klev = fields['KLEV']

        for argname in argnames:
           if argname in argnames_nlev:
               print('Loading field:',argname)
               fields[argname] = np.asfortranarray(np.transpose(np.reshape(
                                   np.ascontiguousarray(f[argname.upper()]),(nblocks,nlev,nproma),order='C')))

           if argname in argnames_nlevp:
               print('Loading field:',argname)
               fields[argname] = np.asfortranarray(np.transpose(np.reshape(
                                   np.ascontiguousarray(f[argname.upper()]),(nblocks,nlev+1,nproma),order='C')))

           if argname in argnames_withnclv:
               print('Loading field:',argname)
               fields[argname] = np.asfortranarray(np.transpose(np.reshape(
                                   np.ascontiguousarray(f[argname.upper()]),(nblocks,NCLV,nlev,nproma),order='C')))

           if argname in argnames_tend:
               print('Loading field:',argname)
               fields[argname] = np.asfortranarray(np.transpose(np.reshape(
                                   np.ascontiguousarray(f[argname.upper()]),(nblocks,nlev,nproma),order='C')))

    argnames_nlev = [
        'pqsat', 'pcovptot'
    ]

    argnames_nlevp = [
        'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn'
    ]

    argnames_buffer = [
        'buffer_cml'
    ]
 
    argnames_tend = [
        'tendency_loc_a','tendency_loc_t','tendency_loc_q',
        'tendency_cml_a'
    ]
    
    argnames_tend_cld = [
        'tendency_loc_cld'
    ]


    for argname in argnames_nlev:
        fields[argname] = np.zeros(shape=(nproma,nlev  ,nblocks), order='F')

    for argname in argnames_nlevp:
        fields[argname] = np.zeros(shape=(nproma,nlev+1,nblocks), order='F')

    for argname in argnames_buffer:
        fields[argname] = np.zeros(shape=(nproma,nlev,3+NCLV,nblocks), order='F')

    for argname in argnames_tend:
        fields[argname] = np.zeros(shape=(nproma,nlev,nblocks), order='F')

    for argname in argnames_tend:
        fields[argname] = np.zeros(shape=(nproma,nlev,NCLV,nblocks), order='F')

    for argname in argnames_buffer:
        fields[argname] = np.zeros(shape=(nproma,nlev,3+NCLV,nblocks), order='F')

    pack_buffer_using_tendencies(fields['buffer_cml'],
                                 fields['tendency_cml_a'],
                                 fields['tendency_cml_t'],
                                 fields['tendency_cml_q'],
                                 fields['tendency_cml_cld'])
    return fields

def pack_buffer_using_tendencies(buffervar,tendency_a,tendency_t,tendency_q,tendency_cld):
     buffervar[:,:,0       ,:]=tendency_t  [:,:,:]
     buffervar[:,:,2       ,:]=tendency_q  [:,:,:]
     buffervar[:,:,3:3+NCLV-1,:]=tendency_cld[:,:,0:NCLV-1,:]

def  unpack_buffer_to_tendencies(buffervar,tendency_a,tendency_t,tendency_q,tendency_cld):
     tendency_t  [:,:,:]=buffervar[:,:,0       ,:]
     tendency_q  [:,:,:]=buffervar[:,:,1       ,:]
     tendency_cld[:,:,0:NCLV-1,:]=buffervar[:,:,3:3+NCLV-1,:]

def load_input_parameters(path,yrecldp,yrephli,yrmcst,yrethf,yrecld):
    with h5py.File(path, 'r') as f:
        tecldp_keys = [k for k in f.keys() if 'YRECLDP' in k]
        for k in tecldp_keys:
            attrkey = k.replace('YRECLDP_', '').lower()
            setattr(yrecldp, attrkey, f[k][0])

        tephli_keys = [k for k in f.keys() if 'YREPHLI' in k]
        for k in tephli_keys:
            attrkey = k.replace('YREPHLI_', '').lower()
            setattr(yrephli, attrkey, f[k][0])

        yrmcst.rg = f['RG'][0]
        yrmcst.rd = f['RD'][0]
        yrmcst.rcpd = f['RCPD'][0]
        yrmcst.retv = f['RETV'][0]
        yrmcst.rlvtt = f['RLVTT'][0]
        yrmcst.rlstt = f['RLSTT'][0]
        yrmcst.rlmlt = f['RLMLT'][0]
        yrmcst.rtt = f['RTT'][0]
        yrmcst.rv = f['RV'][0]

        yrethf.r2es = f['R2ES'][0]
        yrethf.r3les = f['R3LES'][0]
        yrethf.r3ies = f['R3IES'][0]
        yrethf.r4les = f['R4LES'][0]
        yrethf.r4ies = f['R4IES'][0]
        yrethf.r5les = f['R5LES'][0]
        yrethf.r5ies = f['R5IES'][0]
        yrethf.r5alvcp = f['R5ALVCP'][0]
        yrethf.r5alscp = f['R5ALSCP'][0]
        yrethf.ralvdcp = f['RALVDCP'][0]
        yrethf.ralsdcp = f['RALSDCP'][0]
        yrethf.ralfdcp = f['RALFDCP'][0]
        yrethf.rtwat = f['RTWAT'][0]
        yrethf.rtice = f['RTICE'][0]
        yrethf.rticecu = f['RTICECU'][0]
        yrethf.rtwat_rtice_r = f['RTWAT_RTICE_R'][0]
        yrethf.rtwat_rticecu_r = f['RTWAT_RTICECU_R'][0]
        yrethf.rkoop1 = f['RKOOP1'][0]
        yrethf.rkoop2 = f['RKOOP2'][0]

        yrethf.rvtmp2 = 0.0

        klev = f['KLEV'][0]
        pap = np.ascontiguousarray(f['PAP'])
        paph = np.ascontiguousarray(f['PAPH'])
        yrecld.ceta = np.ndarray(order="C", shape=(klev, ))
        yrecld.ceta[:] = pap[0:,0] / paph[klev,0]

        yrephli.lphylin = True

    return yrecldp, yrmcst, yrethf, yrephli, yrecld


def load_reference_fortran_fields(path):
    """

    """
    fields = OrderedDict()

    argnames = [
        'plude', 'pcovptot', 'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn',
        'tendency_loc_a', 'tendency_loc_q', 'tendency_loc_t', 'tendency_loc_cld',
    ]

    with h5py.File(path, 'r') as f:
        for argname in argnames:
            fields[argname.lower()] = np.ascontiguousarray(f[argname])

    return fields



def cloudsc_validate(fields, ref_fields, cloudsc_args):
    # List of refencece fields names in order
    _field_names = [
        'plude', 'pcovptot', 'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn',
        'tendency_loc_a', 'tendency_loc_q', 'tendency_loc_t', # 'tendency_loc_cld',
    ]
    kidia = cloudsc_args['kidia']
    kfdia = cloudsc_args['kfdia']
    ngptot = kfdia - kidia + 1

    print("             Variable Dim             MinValue             MaxValue            AbsMaxErr         AvgAbsErr/GP          MaxRelErr-%")
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
                  maxrel=0.0 if zerrsum < zeps else (zerrsum/(1.0+zsum) if zsum < zeps else zerrsum/zsum)
              )
        )

def convert_fortran_output_to_python (nproma,nlev,nblocks,plude, pcovptot, pfplsl, pfplsn, pfhpsl, pfhpsn, buffer_loc ):
    """

    """

    fields = OrderedDict()
    argnames_nlev = [
        'plude', 'pcovptot'
    ]

    argnames_nlevp = [
        'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn'
    ]

    argnames_buffer = [
        'buffer_loc'
    ]
    argnames_tend = [
        'tendency_loc_a','tendency_loc_t','tendency_loc_q',
    ]
    
    argnames_tend_cld = [
        'tendency_loc_cld'
    ]


    for argname in argnames_nlev:
        fields[argname] = np.ascontiguousarray(np.transpose(locals()[argname][:,:,0]))

    for argname in argnames_nlevp:
        fields[argname] = np.ascontiguousarray(np.transpose(locals()[argname][:,:,0]))

    for argname in argnames_tend:
        locals()[argname] = np.zeros(shape=(nproma,nlev,nblocks), order='F')

    for argname in argnames_tend_cld:
        locals()[argname] = np.zeros(shape=(nproma,nlev,NCLV,nblocks), order='F')


    unpack_buffer_to_tendencies(locals() ['buffer_loc'],
                                locals() ['tendency_loc_a'],
                                locals() ['tendency_loc_t'],
                                locals() ['tendency_loc_q'],
                                locals() ['tendency_loc_cld'])

    fields['tendency_loc_a'] = np.ascontiguousarray(np.transpose(locals()['tendency_loc_a'][:,:,0]))
    fields['tendency_loc_t'] = np.ascontiguousarray(np.transpose(locals()['tendency_loc_t'][:,:,0])) 
    fields['tendency_loc_q'] = np.ascontiguousarray(np.transpose(locals()['tendency_loc_q'][:,:,0])) 
 

    return fields

def load_reference_fields (path):
    """

    """

    fields = OrderedDict()
    argnames_nlev = [
        'plude', 'pcovptot'
    ]

    argnames_nlevp = [
        'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn'
    ]

    argnames_tend = [
        'tendency_loc_a','tendency_loc_t','tendency_loc_q',
    ]
    
    argnames_tend_cld = [
        'tendency_loc_cld'
    ]

    with h5py.File(path, 'r') as f:
        for argname in argnames_nlev:
            fields[argname] = np.ascontiguousarray(f[argname.upper()])
    
        for argname in argnames_nlevp:
            fields[argname] = np.ascontiguousarray(f[argname.upper()])
    
        for argname in argnames_tend:
            fields[argname] = np.ascontiguousarray(f[argname.upper()])
    
        for argname in argnames_tend_cld:
            fields[argname] = np.ascontiguousarray(f[argname.upper()]) 


    return fields

def cloudsc_validate(fields, ref_fields):
    # List of refencece fields names in order
    _field_names = [
        'plude', 'pcovptot', 'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn', 
        'tendency_loc_a', 'tendency_loc_q', 'tendency_loc_t', # 'tendency_loc_cld',
    ]
    kidia = 1 
    kfdia = 100 
    ngptot = kfdia - kidia + 1

    print("             Variable Dim             MinValue             MaxValue            AbsMaxErr         AvgAbsErr/GP          MaxRelErr-%")
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
                  maxrel=0.0 if zerrsum < zeps else (zerrsum/(1.0+zsum) if zsum < zeps else zerrsum/zsum)
              )
        )


