import h5py
import numpy as np
import cloudsc as clsc
from pathlib import Path
from collections import OrderedDict


NCLV = 5      # number of microphysics variables


def define_fortran_fields(nproma,nlev,nblocks):
    """

    """
    fields = OrderedDict()

#   argnames_nlev = [
#       'pt', 'pq',
#       'pvfa', 'pvfl', 'pvfi', 'pdyna', 'pdynl', 'pdyni',
#       'phrsw', 'phrlw','pvervel','pap','plu','plude',
#       'psnde', 'pmfu', 'pmfd',
#       'pa', 'psupsat', 
#       'plcrit_aer','picrit_aer', 'pre_ice', 
#       'pccn', 'pnice',
#       'pcovptot',

#   ]

#   argnames_nlevp = [
#       'paph', 
#       'pfsqlf',   'pfsqif' ,  'pfcqnng',  'pfcqlng',
#       'pfsqrf',   'pfsqsf' ,  'pfcqrng',  'pfcqsng',
#       'pfsqltur', 'pfsqitur' ,
#       'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn'
#   ]

#   argnames_withnclv= [
#       'pclv','tendency_loc_cld'
#   ]

#   argnames_buffer = [
#       'buffer_tmp', 'buffer_loc'
#   ]
#
#   argnames_tend = [
#       'tendency_loc_a','tendency_loc_t','tendency_loc_q'
#   ]
#   
#   argnames_nproma = [
#       'plsm', 'ldcum', 'ktype', 'prainfrac_toprfz'
#   ]

#   for argname in argnames_nlev:
#       fields[argname] = np.zeros(shape=(nproma,nlev  ,nblocks), order='F')

#   for argname in argnames_nlevp:
#       fields[argname] = np.zeros(shape=(nproma,nlev+1,nblocks), order='F')

#   for argname in argnames_withnclv:
#       fields[argname] = np.zeros(shape=(nproma,nlev,NCLV,nblocks), order='F')

#   for argname in argnames_buffer:
#       fields[argname] = np.zeros(shape=(nproma,nlev,3+NCLV,nblocks), order='F')

#   for argname in argnames_tend:
#       fields[argname] = np.zeros(shape=(nproma,nlev,nblocks), order='F')

#   for argname in argnames_tend:
#       fields[argname] = np.zeros(shape=(nproma,nblocks), order='F')

    fields['ydomcst']=clsc.yomcst.TOMCST()
    fields['ydoethf']=clsc.yoethf.TOETHF()
    fields['ydecldp']=clsc.yoecldp.TECLDP()
    fields['ydephli']=clsc.yoephli.TEPHLI()

    return fields

def load_input_fortran_fields(path, nproma, nlev, nblocks,  transpose=False):
    """

    """
    fields = OrderedDict()

    argnames_nlev = [
        'pt', 'pq',
        'pvfa', 'pvfl', 'pvfi', 'pdyna', 'pdynl', 'pdyni',
        'phrsw', 'phrlw','pvervel','pap','plu','plude',
        'psnde', 'pmfu', 'pmfd',
        'pa', 'psupsat', 
        'plcrit_aer','picrit_aer', 'pre_ice', 
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

    argnames = [
        'pt', 'pq',
        'pvfa', 'pvfl', 'pvfi', 'pdyna', 'pdynl', 'pdyni',
        'phrsw', 'phrlw','pvervel','pap','plu','plude',
        'psnde', 'pmfu', 'pmfd',
        'pa', 'psupsat', 
        'plcrit_aer','picrit_aer', 'pre_ice', 
        'pccn', 'pnice',
        'paph',
        'pclv','tendency_tmp_cld',
        'tendency_tmp_t','tendency_tmp_q','tendency_tmp_a',
        'kfldx',
         'plsm', 'ldcum', 'ktype'
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

           if argname in argnames_nproma:
               print('Loading field:',argname)
               fields[argname] = np.asfortranarray(np.transpose(np.reshape(
                                   np.ascontiguousarray(f[argname.upper()]),(nblocks,nproma),order='C')))
           if argname in argnames_scalar:
               print('Loading field:',argname)
               fields[argname] = f[argname.upper()][0]

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
        print('Define null field:',argname)
        fields[argname] = np.zeros(shape=(nproma,nlev  ,nblocks), order='F')

    for argname in argnames_nlevp:
        print('Define null field:',argname)
        fields[argname] = np.zeros(shape=(nproma,nlev+1,nblocks), order='F')

    for argname in argnames_buffer:
        print('Define null field:',argname)
        fields[argname] = np.zeros(shape=(nproma,nlev,3+NCLV,nblocks), order='F')

    for argname in argnames_tend:
        print('Define null field:',argname)
        fields[argname] = np.zeros(shape=(nproma,nlev,nblocks), order='F')

    for argname in argnames_tend:
        print('Define null field:',argname)
        fields[argname] = np.zeros(shape=(nproma,nlev,NCLV,nblocks), order='F')

    for argname in argnames_buffer:
        print('Define null field:',argname)
        fields[argname] = np.zeros(shape=(nproma,nlev,3+NCLV,nblocks), order='F')

    for argname in argnames_nproma:
        print('Define null field:',argname)
        fields[argname] = np.zeros(shape=(nproma,nblocks), order='F')

    pack_buffer_using_tendencies(fields['buffer_tmp'],
                                 fields['tendency_tmp_a'],
                                 fields['tendency_tmp_t'],
                                 fields['tendency_tmp_q'],
                                 fields['tendency_tmp_cld'])
    return fields

def pack_buffer_using_tendencies(buffervar,tendency_a,tendency_t,tendency_q,tendency_cld):
     buffervar[:,:,0       ,:]=tendency_t  [:,:,:]
     buffervar[:,:,1       ,:]=tendency_a  [:,:,:]
     buffervar[:,:,2       ,:]=tendency_q  [:,:,:]
     buffervar[:,:,3:3+NCLV-1,:]=tendency_cld[:,:,0:NCLV-1,:]

def  unpack_buffer_to_tendencies(buffervar,tendency_a,tendency_t,tendency_q,tendency_cld):
     tendency_t  [:,:,:]=buffervar[:,:,0       ,:]
     tendency_a  [:,:,:]=buffervar[:,:,1       ,:]
     tendency_q  [:,:,:]=buffervar[:,:,2       ,:]
     tendency_cld[:,:,0:NCLV-1,:]=buffervar[:,:,3:3+NCLV-1,:]

def load_input_parameters(path,yrecldp,yrephli,yrmcst,yrethf):
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

        yrecldp.laericeauto           = f['YRECLDP_LAERICEAUTO'][0]         
        yrecldp.laericesed            = f['YRECLDP_LAERICESED'][0]
        yrecldp.laerliqautolsp        = f['YRECLDP_LAERLIQAUTOLSP'][0]
        yrecldp.laerliqcoll           = f['YRECLDP_LAERLIQCOLL'][0]
        yrecldp.lcldbudget            = f['YRECLDP_LCLDBUDGET'][0]
        yrecldp.ncldtop               = f['YRECLDP_NCLDTOP'][0]
        yrecldp.nssopt                = f['YRECLDP_NSSOPT'][0]
        yrecldp.ramid                 = f['YRECLDP_RAMID'][0]
        yrecldp.ramin                 = f['YRECLDP_RAMIN'][0]
        yrecldp.rccn                  = f['YRECLDP_RCCN'][0]
        yrecldp.rclcrit_land          = f['YRECLDP_RCLCRIT_LAND'][0]
        yrecldp.rclcrit_sea           = f['YRECLDP_RCLCRIT_SEA'][0]
        yrecldp.rcldiff               = f['YRECLDP_RCLDIFF'][0]
        yrecldp.rcldiff_convi         = f['YRECLDP_RCLDIFF_CONVI'][0]
        yrecldp.rcldtopcf             = f['YRECLDP_RCLDTOPCF'][0]
        yrecldp.rcl_apb1              = f['YRECLDP_RCL_APB1'][0]
        yrecldp.rcl_apb2              = f['YRECLDP_RCL_APB2'][0]
        yrecldp.rcl_apb3              = f['YRECLDP_RCL_APB3'][0]
        yrecldp.rcl_cdenom1           = f['YRECLDP_RCL_CDENOM1'][0]
        yrecldp.rcl_cdenom2           = f['YRECLDP_RCL_CDENOM2'][0]
        yrecldp.rcl_cdenom3           = f['YRECLDP_RCL_CDENOM3'][0]
        yrecldp.rcl_const1i           = f['YRECLDP_RCL_CONST1I'][0]
        yrecldp.rcl_const1r           = f['YRECLDP_RCL_CONST1R'][0]
        yrecldp.rcl_const1s           = f['YRECLDP_RCL_CONST1S'][0]
        yrecldp.rcl_const2i           = f['YRECLDP_RCL_CONST2I'][0]
        yrecldp.rcl_const2r           = f['YRECLDP_RCL_CONST2R'][0]
        yrecldp.rcl_const2s           = f['YRECLDP_RCL_CONST2S'][0]
        yrecldp.rcl_const3i           = f['YRECLDP_RCL_CONST3I'][0]
        yrecldp.rcl_const3r           = f['YRECLDP_RCL_CONST3R'][0]
        yrecldp.rcl_const3s           = f['YRECLDP_RCL_CONST3S'][0]
        yrecldp.rcl_const4i           = f['YRECLDP_RCL_CONST4I'][0]
        yrecldp.rcl_const4r           = f['YRECLDP_RCL_CONST4R'][0]
        yrecldp.rcl_const4s           = f['YRECLDP_RCL_CONST4S'][0]
        yrecldp.rcl_const5i           = f['YRECLDP_RCL_CONST5I'][0]
        yrecldp.rcl_const5r           = f['YRECLDP_RCL_CONST5R'][0]
        yrecldp.rcl_const5s           = f['YRECLDP_RCL_CONST5S'][0]
        yrecldp.rcl_const6i           = f['YRECLDP_RCL_CONST6I'][0]
        yrecldp.rcl_const6r           = f['YRECLDP_RCL_CONST6R'][0]
        yrecldp.rcl_const6s           = f['YRECLDP_RCL_CONST6S'][0]
        yrecldp.rcl_const7s           = f['YRECLDP_RCL_CONST7S'][0]
        yrecldp.rcl_const8s           = f['YRECLDP_RCL_CONST8S'][0]
        yrecldp.rcl_fac1              = f['YRECLDP_RCL_FAC1'][0]
        yrecldp.rcl_fac2              = f['YRECLDP_RCL_FAC2'][0]
        yrecldp.rcl_fzrab             = f['YRECLDP_RCL_FZRAB'][0]
        yrecldp.rcl_ka273             = f['YRECLDP_RCL_KA273'][0]
        yrecldp.rcl_kkaac             = f['YRECLDP_RCL_KKAac'][0]
        yrecldp.rcl_kkaau             = f['YRECLDP_RCL_KKAau'][0]
        yrecldp.rcl_kkbac             = f['YRECLDP_RCL_KKBac'][0]
        yrecldp.rcl_kkbaun            = f['YRECLDP_RCL_KKBaun'][0]
        yrecldp.rcl_kkbauq            = f['YRECLDP_RCL_KKBauq'][0]
        yrecldp.rcl_kk_cloud_num_land = f['YRECLDP_RCL_KK_cloud_num_land'][0]
        yrecldp.rcl_kk_cloud_num_sea  = f['YRECLDP_RCL_KK_cloud_num_sea'][0]
        yrecldp.rcl_x3i               = f['YRECLDP_RCL_X3I'][0]
        yrecldp.rcovpmin              = f['YRECLDP_RCOVPMIN'][0]
        yrecldp.rdensref              = f['YRECLDP_RDENSREF'][0]
        yrecldp.rdepliqrefdepth       = f['YRECLDP_RDEPLIQREFDEPTH'][0]
        yrecldp.rdepliqrefrate        = f['YRECLDP_RDEPLIQREFRATE'][0]
        yrecldp.ricehi1               = f['YRECLDP_RICEHI1'][0]
        yrecldp.ricehi2               = f['YRECLDP_RICEHI2'][0]
        yrecldp.riceinit              = f['YRECLDP_RICEINIT'][0]
        yrecldp.rkconv                = f['YRECLDP_RKCONV'][0]
        yrecldp.rkooptau              = f['YRECLDP_RKOOPTAU'][0]
        yrecldp.rlcritsnow            = f['YRECLDP_RLCRITSNOW'][0]
        yrecldp.rlmin                 = f['YRECLDP_RLMIN'][0]
        yrecldp.rnice                 = f['YRECLDP_RNICE'][0]
        yrecldp.rpecons               = f['YRECLDP_RPECONS'][0]
        yrecldp.rprc1                 = f['YRECLDP_RPRC1'][0]
        yrecldp.rprecrhmax            = f['YRECLDP_RPRECRHMAX'][0]
        yrecldp.rsnowlin1             = f['YRECLDP_RSNOWLIN1'][0]
        yrecldp.rsnowlin2             = f['YRECLDP_RSNOWLIN2'][0]
        yrecldp.rtaumel               = f['YRECLDP_RTAUMEL'][0]
        yrecldp.rthomo                = f['YRECLDP_RTHOMO'][0]
        yrecldp.rvice                 = f['YRECLDP_RVICE'][0]
        yrecldp.rvrain                = f['YRECLDP_RVRAIN'][0]
        yrecldp.rvrfactor             = f['YRECLDP_RVRFACTOR'][0]
        yrecldp.rvsnow                = f['YRECLDP_RVSNOW'][0]
        klev = f['KLEV'][0]
        pap = np.ascontiguousarray(f['PAP'])
        paph = np.ascontiguousarray(f['PAPH'])

        yrephli.lphylin = True

    return yrecldp, yrmcst, yrethf, yrephli


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

def convert_fortran_output_to_python (nproma,nlev,nblocks,
                                      plude, pcovptot, 
                                      pfplsl, pfplsn, pfhpsl, pfhpsn,
                                      pfsqlf, pfsqif, pfcqnng,  pfcqlng,
                                      pfsqrf,   pfsqsf ,  pfcqrng,  pfcqsng,
                                      pfsqltur, pfsqitur, 
                                      prainfrac_toprfz,
                                      buffer_loc ):
    """

    """

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

    argnames_buffer = [
        'buffer_loc'
    ]
    argnames_tend = [
        'tendency_loc_a','tendency_loc_t','tendency_loc_q'
    ]
    
    argnames_tend_cld = [
        'tendency_loc_cld'
    ]
    
    argnames_nproma = [ 
        'prainfrac_toprfz'
    ] 

    for argname in argnames_nlev:
        fields[argname] = np.ascontiguousarray(np.transpose(locals()[argname][:,:,0]))

    for argname in argnames_nlevp:
        fields[argname] = np.ascontiguousarray(np.transpose(locals()[argname][:,:,0]))

    for argname in argnames_nproma:
        fields[argname] = np.ascontiguousarray(np.transpose(locals()[argname][:,0]))

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
        for argname in argnames_nlev:
            fields[argname] = np.ascontiguousarray(f[argname.upper()])
    
        for argname in argnames_nlevp:
            fields[argname] = np.ascontiguousarray(f[argname.upper()])
    
        for argname in argnames_nproma:
            fields[argname] = np.ascontiguousarray(f[argname.upper()])

        for argname in argnames_tend:
            fields[argname] = np.ascontiguousarray(f[argname.upper()])
    
        for argname in argnames_tend_cld:
            fields[argname] = np.ascontiguousarray(f[argname.upper()]) 


    return fields

def cloudsc_validate(fields, ref_fields):
    # List of refencece fields names in order
    _field_names = [
        'plude', 'pcovptot','prainfrac_toprfz',
        'pfsqlf',   'pfsqif' ,  'pfcqnng',  'pfcqlng',
        'pfsqrf',   'pfsqsf' ,  'pfcqrng',  'pfcqsng',
        'pfsqltur', 'pfsqitur' ,
        'pfplsl', 'pfplsn', 'pfhpsl', 'pfhpsn', 
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


