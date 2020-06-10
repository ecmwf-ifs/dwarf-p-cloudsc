MODULE cloudsc_claw_mod

CONTAINS
 SUBROUTINE cloudsc_claw ( klev , ptsphy , pt , pq , tendency_cml_t , tendency_cml_q , tendency_cml_cld , tendency_tmp_t ,&
  tendency_tmp_q , tendency_tmp_a , tendency_tmp_cld , tendency_loc_t , tendency_loc_q , tendency_loc_a , tendency_loc_cld , pvfa&
  , pvfl , pvfi , pdyna , pdynl , pdyni , phrsw , phrlw , pvervel , pap , paph , plsm , ldcum , ktype , plu , plude , psnde , pmfu&
  , pmfd , pa , pclv , psupsat , plcrit_aer , picrit_aer , pre_ice , pccn , pnice , pcovptot , prainfrac_toprfz , pfsqlf , pfsqif&
  , pfcqnng , pfcqlng , pfsqrf , pfsqsf , pfcqrng , pfcqsng , pfsqltur , pfsqitur , pfplsl , pfplsn , pfhpsl , pfhpsn ,&
  kfldx , yrecldp , nproma )
  USE parkind1 , ONLY: jpim , jprb
  USE yomcst , ONLY: rg , rd , rcpd , retv , rlvtt , rlstt , rlmlt , rtt , rv
  USE yoethf , ONLY: r2es , r3les , r3ies , r4les , r4ies , r5les , r5ies , r5alvcp , r5alscp , ralvdcp , ralsdcp , ralfdcp ,&
   rtwat , rtice , rticecu , rtwat_rtice_r , rtwat_rticecu_r , rkoop1 , rkoop2
  USE yoecldp , ONLY: tecldp , ncldqv , ncldql , ncldqr , ncldqi , ncldqs , nclv
  USE yomphyder , ONLY: state_type
  INTEGER , INTENT(IN) :: nproma
  INTERFACE
   SUBROUTINE abor1 ( cdtext )

    CHARACTER ( LEN= * ) :: cdtext
   END SUBROUTINE abor1
  END INTERFACE

  REAL ( KIND=JPRB ) , INTENT(IN) :: plcrit_aer ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: picrit_aer ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: pre_ice ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: pccn ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: pnice ( 1 : nproma , 1 : klev )
  INTEGER ( KIND= 4 ) , INTENT(IN) :: klev
  REAL ( KIND=JPRB ) , INTENT(IN) :: ptsphy
  REAL ( KIND=JPRB ) , INTENT(IN) :: pt ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: pq ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: tendency_cml_t ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: tendency_cml_q ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: tendency_cml_cld ( 1 : nproma , 1 : klev , 1 : 5 )
  REAL ( KIND=JPRB ) , INTENT(IN) :: tendency_tmp_t ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: tendency_tmp_q ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: tendency_tmp_a ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: tendency_tmp_cld ( 1 : nproma , 1 : klev , 1 : 5 )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: tendency_loc_t ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: tendency_loc_q ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: tendency_loc_a ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: tendency_loc_cld ( 1 : nproma , 1 : klev , 1 : 5 )
  REAL ( KIND=JPRB ) , INTENT(IN) :: pvfa ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: pvfl ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: pvfi ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: pdyna ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: pdynl ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: pdyni ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: phrsw ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: phrlw ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: pvervel ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: pap ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: paph ( 1 : nproma , 1 : klev + 1 )
  REAL ( KIND=JPRB ) , INTENT(IN) :: plsm ( 1 : nproma )
  LOGICAL , INTENT(IN) :: ldcum ( 1 : nproma )
  INTEGER ( KIND= 4 ) , INTENT(IN) :: ktype ( 1 : nproma )
  REAL ( KIND=JPRB ) , INTENT(IN) :: plu ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(INOUT) :: plude ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: psnde ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: pmfu ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: pmfd ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(IN) :: pa ( 1 : nproma , 1 : klev )
  INTEGER ( KIND= 4 ) , INTENT(IN) :: kfldx
  REAL ( KIND=JPRB ) , INTENT(IN) :: pclv ( 1 : nproma , 1 : klev , 1 : 5 )
  REAL ( KIND=JPRB ) , INTENT(IN) :: psupsat ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: pcovptot ( 1 : nproma , 1 : klev )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: prainfrac_toprfz ( 1 : nproma )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: pfsqlf ( 1 : nproma , 1 : klev + 1 )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: pfsqif ( 1 : nproma , 1 : klev + 1 )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: pfcqlng ( 1 : nproma , 1 : klev + 1 )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: pfcqnng ( 1 : nproma , 1 : klev + 1 )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: pfsqrf ( 1 : nproma , 1 : klev + 1 )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: pfsqsf ( 1 : nproma , 1 : klev + 1 )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: pfcqrng ( 1 : nproma , 1 : klev + 1 )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: pfcqsng ( 1 : nproma , 1 : klev + 1 )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: pfsqltur ( 1 : nproma , 1 : klev + 1 )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: pfsqitur ( 1 : nproma , 1 : klev + 1 )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: pfplsl ( 1 : nproma , 1 : klev + 1 )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: pfplsn ( 1 : nproma , 1 : klev + 1 )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: pfhpsl ( 1 : nproma , 1 : klev + 1 )
  REAL ( KIND=JPRB ) , INTENT(OUT) :: pfhpsn ( 1 : nproma , 1 : klev + 1 )
  TYPE ( tecldp ) , INTENT(INOUT) :: yrecldp
  REAL ( KIND=JPRB ) :: zlcond1 ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zlcond2 ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zlevap
  REAL ( KIND=JPRB ) :: zleros
  REAL ( KIND=JPRB ) :: zlevapl ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zlevapi ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zrainaut ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zsnowaut ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zliqcld ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zicecld ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zfokoop ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zfoealfa ( 1 : klev + 1 )
  REAL ( KIND=JPRB ) :: zicenuclei ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zlicld ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zacond
  REAL ( KIND=JPRB ) :: zaeros
  REAL ( KIND=JPRB ) :: zlfinalsum ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zdqs ( 1 : nproma )
  REAL ( KIND=JPRB ) :: ztold ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zqold ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zdtgdp ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zrdtgdp ( 1 : nproma )
  REAL ( KIND=JPRB ) :: ztrpaus ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zcovpclr ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zpreclr
  REAL ( KIND=JPRB ) :: zcovptot ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zcovpmax ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zqpretot ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zdpevap
  REAL ( KIND=JPRB ) :: zdtforc
  REAL ( KIND=JPRB ) :: zdtdiab
  REAL ( KIND=JPRB ) :: ztp1 ( 1 : klev )
  REAL ( KIND=JPRB ) :: zldefr ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zldifdt ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zdtgdpf ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zlcust ( 1 : 5 )
  REAL ( KIND=JPRB ) :: zacust ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zmf ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zrho ( 1 : nproma )
  REAL ( KIND=JPRB ) :: ztmp1 ( 1 : nproma )
  REAL ( KIND=JPRB ) :: ztmp2 ( 1 : nproma )
  REAL ( KIND=JPRB ) :: ztmp3 ( 1 : nproma )
  REAL ( KIND=JPRB ) :: ztmp4 ( 1 : nproma )
  REAL ( KIND=JPRB ) :: ztmp5 ( 1 : nproma )
  REAL ( KIND=JPRB ) :: ztmp6 ( 1 : nproma )
  REAL ( KIND=JPRB ) :: ztmp7 ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zalfawm ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zsolab ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zsolac ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zanew
  REAL ( KIND=JPRB ) :: zanewm1 ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zgdp ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zda ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zli ( 1 : klev )
  REAL ( KIND=JPRB ) :: za ( 1 : klev )
  REAL ( KIND=JPRB ) :: zaorig ( 1 : klev )
  LOGICAL :: llflag
  LOGICAL :: llo1
  INTEGER ( KIND= 4 ) :: icall
  INTEGER ( KIND= 4 ) :: ik
  INTEGER ( KIND= 4 ) :: jk
  INTEGER ( KIND= 4 ) :: jm
  INTEGER ( KIND= 4 ) :: jn
  INTEGER ( KIND= 4 ) :: jo
  INTEGER ( KIND= 4 ) :: jlen
  INTEGER ( KIND= 4 ) :: is
  REAL ( KIND=JPRB ) :: zdp ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zpaphd ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zalfa
  REAL ( KIND=JPRB ) :: zalfaw
  REAL ( KIND=JPRB ) :: zbeta
  REAL ( KIND=JPRB ) :: zbeta1
  REAL ( KIND=JPRB ) :: zcfpr
  REAL ( KIND=JPRB ) :: zcor
  REAL ( KIND=JPRB ) :: zcdmax
  REAL ( KIND=JPRB ) :: zmin ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zlcondlim
  REAL ( KIND=JPRB ) :: zdenom
  REAL ( KIND=JPRB ) :: zdpmxdt
  REAL ( KIND=JPRB ) :: zdpr
  REAL ( KIND=JPRB ) :: zdtdp
  REAL ( KIND=JPRB ) :: ze
  REAL ( KIND=JPRB ) :: zepsec
  REAL ( KIND=JPRB ) :: zfac
  REAL ( KIND=JPRB ) :: zfaci
  REAL ( KIND=JPRB ) :: zfacw
  REAL ( KIND=JPRB ) :: zgdcp
  REAL ( KIND=JPRB ) :: zinew
  REAL ( KIND=JPRB ) :: zlcrit
  REAL ( KIND=JPRB ) :: zmfdn
  REAL ( KIND=JPRB ) :: zprecip
  REAL ( KIND=JPRB ) :: zqe
  REAL ( KIND=JPRB ) :: zqsat
  REAL ( KIND=JPRB ) :: zqtmst
  REAL ( KIND=JPRB ) :: zrdcp
  REAL ( KIND=JPRB ) :: zrhc
  REAL ( KIND=JPRB ) :: zsig
  REAL ( KIND=JPRB ) :: zsigk
  REAL ( KIND=JPRB ) :: zwtot
  REAL ( KIND=JPRB ) :: zzco
  REAL ( KIND=JPRB ) :: zzdl
  REAL ( KIND=JPRB ) :: zzrh
  REAL ( KIND=JPRB ) :: zzzdt
  REAL ( KIND=JPRB ) :: zqadj
  REAL ( KIND=JPRB ) :: zqnew
  REAL ( KIND=JPRB ) :: ztnew
  REAL ( KIND=JPRB ) :: zrg_r
  REAL ( KIND=JPRB ) :: zgdph_r
  REAL ( KIND=JPRB ) :: zcons1
  REAL ( KIND=JPRB ) :: zcond
  REAL ( KIND=JPRB ) :: zcons1a
  REAL ( KIND=JPRB ) :: zlfinal
  REAL ( KIND=JPRB ) :: zmelt
  REAL ( KIND=JPRB ) :: zevap
  REAL ( KIND=JPRB ) :: zfrz
  REAL ( KIND=JPRB ) :: zvpliq
  REAL ( KIND=JPRB ) :: zvpice
  REAL ( KIND=JPRB ) :: zadd
  REAL ( KIND=JPRB ) :: zbdd
  REAL ( KIND=JPRB ) :: zcvds
  REAL ( KIND=JPRB ) :: zice0
  REAL ( KIND=JPRB ) :: zdepos
  REAL ( KIND=JPRB ) :: zsupsat ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zfall
  REAL ( KIND=JPRB ) :: zre_ice
  REAL ( KIND=JPRB ) :: zrldcp
  REAL ( KIND=JPRB ) :: zqp1env
  INTEGER ( KIND= 4 ) :: iphase ( 1 : 5 )
  INTEGER ( KIND= 4 ) :: imelt ( 1 : 5 )
  LOGICAL :: llfall ( 1 : 5 )
  LOGICAL :: llindex1 ( 1 : 5 )
  LOGICAL :: llindex3 ( 1 : 5 , 1 : 5 )
  REAL ( KIND=JPRB ) :: zmax
  REAL ( KIND=JPRB ) :: zrat
  INTEGER ( KIND= 4 ) :: iorder ( 1 : 5 )
  REAL ( KIND=JPRB ) :: zliqfrac ( 1 : klev )
  REAL ( KIND=JPRB ) :: zicefrac ( 1 : klev )
  REAL ( KIND=JPRB ) :: zqx ( 1 : klev , 1 : 5 )
  REAL ( KIND=JPRB ) :: zqx0 ( 1 : klev , 1 : 5 )
  REAL ( KIND=JPRB ) :: zqxn ( 1 : 5 )
  REAL ( KIND=JPRB ) :: zqxfg ( 1 : 5 )
  REAL ( KIND=JPRB ) :: zqxnm1 ( 1 : 5 )
  REAL ( KIND=JPRB ) :: zfluxq ( 1 : 5 )
  REAL ( KIND=JPRB ) :: zpfplsx ( 1 : klev + 1 , 1 : 5 )
  REAL ( KIND=JPRB ) :: zlneg ( 1 : klev , 1 : 5 )
  REAL ( KIND=JPRB ) :: zmeltmax ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zfrzmax ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zicetot ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zqxn2d ( 1 : klev , 1 : 5 )
  REAL ( KIND=JPRB ) :: zqsmix ( 1 : klev )
  REAL ( KIND=JPRB ) :: zqsliq ( 1 : klev )
  REAL ( KIND=JPRB ) :: zqsice ( 1 : klev )
  REAL ( KIND=JPRB ) :: zfoeewmt ( 1 : klev )
  REAL ( KIND=JPRB ) :: zfoeew ( 1 : klev )
  REAL ( KIND=JPRB ) :: zfoeeliqt ( 1 : klev )
  REAL ( KIND=JPRB ) :: zdqsliqdt ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zdqsicedt ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zdqsmixdt ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zcorqsliq ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zcorqsice ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zcorqsmix ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zevaplimliq ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zevaplimice ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zevaplimmix ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zsolqa ( 1 : 5 , 1 : 5 )
  REAL ( KIND=JPRB ) :: zsolqb ( 1 : 5 , 1 : 5 )
  REAL ( KIND=JPRB ) :: zqlhs ( 1 : 5 , 1 : 5 )
  REAL ( KIND=JPRB ) :: zvqx ( 1 : 5 )
  REAL ( KIND=JPRB ) :: zexplicit
  REAL ( KIND=JPRB ) :: zratio ( 1 : 5 )
  REAL ( KIND=JPRB ) :: zsinksum ( 1 : 5 )
  REAL ( KIND=JPRB ) :: zfallsink ( 1 : 5 )
  REAL ( KIND=JPRB ) :: zfallsrce ( 1 : 5 )
  REAL ( KIND=JPRB ) :: zconvsrce ( 1 : 5 )
  REAL ( KIND=JPRB ) :: zconvsink ( 1 : 5 )
  REAL ( KIND=JPRB ) :: zpsupsatsrce ( 1 : 5 )
  REAL ( KIND=JPRB ) , PARAMETER :: ztw1 = 1329.31000000000_jprb
  REAL ( KIND=JPRB ) , PARAMETER :: ztw2 = 0.00746150000000000_jprb
  REAL ( KIND=JPRB ) , PARAMETER :: ztw3 = 85000.0000000000_jprb
  REAL ( KIND=JPRB ) , PARAMETER :: ztw4 = 40.6370000000000_jprb
  REAL ( KIND=JPRB ) , PARAMETER :: ztw5 = 275.000000000000_jprb
  REAL ( KIND=JPRB ) :: zsubsat
  REAL ( KIND=JPRB ) :: ztdmtw0
  REAL ( KIND=JPRB ) :: ztcg
  REAL ( KIND=JPRB ) :: zfacx1i
  REAL ( KIND=JPRB ) :: zfacx1s
  REAL ( KIND=JPRB ) :: zaplusb
  REAL ( KIND=JPRB ) :: zcorrfac
  REAL ( KIND=JPRB ) :: zcorrfac2
  REAL ( KIND=JPRB ) :: zpr02
  REAL ( KIND=JPRB ) :: zterm1
  REAL ( KIND=JPRB ) :: zterm2
  REAL ( KIND=JPRB ) :: zcldtopdist ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zinfactor
  INTEGER ( KIND= 4 ) :: iwarmrain
  INTEGER ( KIND= 4 ) :: ievaprain
  INTEGER ( KIND= 4 ) :: ievapsnow
  INTEGER ( KIND= 4 ) :: idepice
  REAL ( KIND=JPRB ) :: zrainacc ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zraincld ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zsnowrime ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zsnowcld ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zesatliq
  REAL ( KIND=JPRB ) :: zfallcorr
  REAL ( KIND=JPRB ) :: zlambda
  REAL ( KIND=JPRB ) :: zevap_denom
  REAL ( KIND=JPRB ) :: zcorr2
  REAL ( KIND=JPRB ) :: zka
  REAL ( KIND=JPRB ) :: zconst
  REAL ( KIND=JPRB ) :: ztemp
  REAL ( KIND=JPRB ) :: zsumq0 ( 1 : klev )
  REAL ( KIND=JPRB ) :: zsumq1 ( 1 : klev )
  REAL ( KIND=JPRB ) :: zerrorq ( 1 : klev )
  REAL ( KIND=JPRB ) :: zsumh0 ( 1 : klev )
  REAL ( KIND=JPRB ) :: zsumh1 ( 1 : klev )
  REAL ( KIND=JPRB ) :: zerrorh ( 1 : klev )
  REAL ( KIND=JPRB ) :: zrain
  REAL ( KIND=JPRB ) :: z_tmp1 ( 1 : nproma )
  REAL ( KIND=JPRB ) :: z_tmp2 ( 1 : nproma )
  REAL ( KIND=JPRB ) :: z_tmp3 ( 1 : nproma )
  REAL ( KIND=JPRB ) :: z_tmp4 ( 1 : nproma )
  REAL ( KIND=JPRB ) :: z_tmp6 ( 1 : nproma )
  REAL ( KIND=JPRB ) :: z_tmp7 ( 1 : nproma )
  REAL ( KIND=JPRB ) :: z_tmpk ( 1 : klev )
  REAL ( KIND=JPRB ) :: zhook_handle
  REAL ( KIND=JPRB ) :: ztmpl
  REAL ( KIND=JPRB ) :: ztmpi
  REAL ( KIND=JPRB ) :: ztmpa
  REAL ( KIND=JPRB ) :: zmm
  REAL ( KIND=JPRB ) :: zrr
  REAL ( KIND=JPRB ) :: zrg ( 1 : nproma )
  REAL ( KIND=JPRB ) :: zbudcc ( 1 : kfldx )
  REAL ( KIND=JPRB ) :: zbudl ( 1 : kfldx )
  REAL ( KIND=JPRB ) :: zbudi ( 1 : kfldx )
  REAL ( KIND=JPRB ) :: zzsum
  REAL ( KIND=JPRB ) :: zzratio
  REAL ( KIND=JPRB ) :: zepsilon
  REAL ( KIND=JPRB ) :: zcond1
  REAL ( KIND=JPRB ) :: zqp
  REAL ( KIND=JPRB ) :: psum_solqa ( 1 : nproma )
  INTEGER :: jl, j, i

!$acc data &
!$acc present(plcrit_aer,picrit_aer,pre_ice,pccn,pnice,pt,pq,tendency_cml_t,tendency_cml_q,tendency_cml_cld,tendency_tmp_t &
!$acc ,tendency_tmp_q,tendency_tmp_a,tendency_tmp_cld,tendency_loc_t,tendency_loc_q,tendency_loc_a,tendency_loc_cld,pvfa,pvfl &
!$acc ,pvfi,pdyna,pdynl,pdyni,phrsw,phrlw,pvervel,pap,paph,plsm,ldcum,ktype,plu,plude,psnde,pmfu,pmfd,pa,pclv,psupsat &
!$acc ,pcovptot,prainfrac_toprfz,pfsqlf,pfsqif,pfcqlng,pfcqnng,pfsqrf,pfsqsf,pfcqrng,pfcqsng,pfsqltur,pfsqitur,pfplsl,pfplsn &
!$acc ,pfhpsl,pfhpsn)

!$acc parallel
!$acc loop gang

  DO i = 1 , nproma / 64 

!$acc loop worker vector &
!$acc private(zfoealfa,ztp1,zlcust,zli,za,zaorig,iphase,imelt,llfall,llindex1,llindex3,iorder,zliqfrac,zicefrac,zqx,zqx0,zqxn &
!$acc ,zqxfg,zqxnm1,zfluxq,zpfplsx,zlneg,zqxn2d,zqsmix,zqsliq,zqsice,zfoeewmt,zfoeew,zfoeeliqt,zsolqa,zsolqb,zqlhs,zvqx,zratio &
!$acc ,zsinksum,zfallsink,zfallsrce,zconvsrce,zconvsink,zpsupsatsrce,zsumq0,zsumq1,zerrorq,zsumh0,zsumh1,zerrorh,z_tmpk,zbudcc &
!$acc ,zbudl,zbudi)

  DO jl = (i-1)*64 +1 , i*64

   ASSOCIATE ( laericeauto => yrecldp % laericeauto , laericesed => yrecldp % laericesed , laerliqautolsp => yrecldp %&
     laerliqautolsp , laerliqcoll => yrecldp % laerliqcoll , lcldbudget => yrecldp % lcldbudget , ncldtop => yrecldp % ncldtop, &
     nssopt => yrecldp % nssopt , ramid => yrecldp % ramid , ramin => yrecldp % ramin , rccn => yrecldp % rccn , rclcrit_land =>&
     yrecldp % rclcrit_land , rclcrit_sea => yrecldp % rclcrit_sea , rcldiff => yrecldp % rcldiff , rcldiff_convi => yrecldp %&
     rcldiff_convi , rcldtopcf => yrecldp % rcldtopcf , rcl_apb1 => yrecldp % rcl_apb1 , rcl_apb2 => yrecldp % rcl_apb2 , rcl_apb3&
     => yrecldp % rcl_apb3 , rcl_cdenom1 => yrecldp % rcl_cdenom1 , rcl_cdenom2 => yrecldp % rcl_cdenom2 , rcl_cdenom3 => yrecldp&
     % rcl_cdenom3 , rcl_const1i => yrecldp % rcl_const1i , rcl_const1r => yrecldp % rcl_const1r , rcl_const1s => yrecldp %&
     rcl_const1s , rcl_const2i => yrecldp % rcl_const2i , rcl_const2r => yrecldp % rcl_const2r , rcl_const2s => yrecldp %&
     rcl_const2s , rcl_const3i => yrecldp % rcl_const3i , rcl_const3r => yrecldp % rcl_const3r , rcl_const3s => yrecldp %&
     rcl_const3s , rcl_const4i => yrecldp % rcl_const4i , rcl_const4r => yrecldp % rcl_const4r , rcl_const4s => yrecldp %&
     rcl_const4s , rcl_const5i => yrecldp % rcl_const5i , rcl_const5r => yrecldp % rcl_const5r , rcl_const5s => yrecldp %&
     rcl_const5s , rcl_const6i => yrecldp % rcl_const6i , rcl_const6r => yrecldp % rcl_const6r , rcl_const6s => yrecldp %&
     rcl_const6s , rcl_const7s => yrecldp % rcl_const7s , rcl_const8s => yrecldp % rcl_const8s , rcl_fac1 => yrecldp % rcl_fac1 ,&
     rcl_fac2 => yrecldp % rcl_fac2 , rcl_fzrab => yrecldp % rcl_fzrab , rcl_ka273 => yrecldp % rcl_ka273 , rcl_kkaac => yrecldp %&
     rcl_kkaac , rcl_kkaau => yrecldp % rcl_kkaau , rcl_kkbac => yrecldp % rcl_kkbac , rcl_kkbaun => yrecldp % rcl_kkbaun ,&
     rcl_kkbauq => yrecldp % rcl_kkbauq , rcl_kk_cloud_num_land => yrecldp % rcl_kk_cloud_num_land , rcl_kk_cloud_num_sea =>&
     yrecldp % rcl_kk_cloud_num_sea , rcl_x3i => yrecldp % rcl_x3i , rcovpmin => yrecldp % rcovpmin , rdensref => yrecldp %&
     rdensref , rdepliqrefdepth => yrecldp % rdepliqrefdepth , rdepliqrefrate => yrecldp % rdepliqrefrate , ricehi1 => yrecldp %&
     ricehi1 , ricehi2 => yrecldp % ricehi2 , riceinit => yrecldp % riceinit , rkconv => yrecldp % rkconv , rkooptau => yrecldp %&
     rkooptau , rlcritsnow => yrecldp % rlcritsnow , rlmin => yrecldp % rlmin , rnice => yrecldp % rnice , rpecons => yrecldp %&
     rpecons , rprc1 => yrecldp % rprc1 , rprecrhmax => yrecldp % rprecrhmax , rsnowlin1 => yrecldp % rsnowlin1 , rsnowlin2 =>&
     yrecldp % rsnowlin2 , rtaumel => yrecldp % rtaumel , rthomo => yrecldp % rthomo , rvice => yrecldp % rvice , rvrain =>&
     yrecldp % rvrain , rvrfactor => yrecldp % rvrfactor , rvsnow => yrecldp % rvsnow )

    zepsilon = 100.0_jprb * epsilon ( zepsilon )
    llcldbudcc = .FALSE.
    llcldbudl = .FALSE.
    llcldbudi = .FALSE.
    iwarmrain = 2
    ievaprain = 2
    ievapsnow = 1
    idepice = 1
    zqtmst = 1.0_jprb / ptsphy
    zgdcp = rg / rcpd
    zrdcp = rd / rcpd
    zcons1a = rcpd / ( rlmlt * rg * rtaumel )
    zepsec = 1.0e-14_jprb
    zrg_r = 1.0_jprb / rg
    zrldcp = 1.0_jprb / ( ralsdcp - ralvdcp )
    iphase ( 5 ) = 0
    iphase ( 1 ) = 1
    iphase ( 3 ) = 1
    iphase ( 2 ) = 2
    iphase ( 4 ) = 2
    imelt ( 5 ) = ( - 99 )
    imelt ( 1 ) = 2
    imelt ( 3 ) = 4
    imelt ( 2 ) = 3
    imelt ( 4 ) = 3
!$acc loop seq
    DO jk = 1 , klev , 1
     tendency_loc_t ( jl , jk ) = 0.0_jprb
     tendency_loc_q ( jl , jk ) = 0.0_jprb
     tendency_loc_a ( jl , jk ) = 0.0_jprb
    END DO
!$acc loop seq
    DO jm = 1 , 4 , 1
!$acc loop seq
     DO jk = 1 , klev , 1
      tendency_loc_cld ( jl , jk , jm ) = 0.0_jprb
     END DO
    END DO
    zvqx ( 5 ) = 0.0_jprb
    zvqx ( 1 ) = 0.0_jprb
    zvqx ( 2 ) = rvice
    zvqx ( 3 ) = rvrain
    zvqx ( 4 ) = rvsnow
    llfall ( : ) = .FALSE.
!$acc loop seq
    DO jm = 1 , 5 , 1
     IF ( zvqx ( jm ) > 0.0_jprb ) THEN
      llfall ( jm ) = .TRUE.
     END IF
    END DO
    llfall ( 2 ) = .FALSE.
!$acc loop seq
    DO jk = 1 , klev , 1
     ztp1 ( jk ) = pt ( jl , jk ) + ptsphy * tendency_tmp_t ( jl , jk )
     zqx ( jk , 5 ) = pq ( jl , jk ) + ptsphy * tendency_tmp_q ( jl , jk )
     zqx0 ( jk , 5 ) = pq ( jl , jk ) + ptsphy * tendency_tmp_q ( jl , jk )
     za ( jk ) = pa ( jl , jk ) + ptsphy * tendency_tmp_a ( jl , jk )
     zaorig ( jk ) = pa ( jl , jk ) + ptsphy * tendency_tmp_a ( jl , jk )
    END DO
!$acc loop seq
    DO jm = 1 , 4 , 1
!$acc loop seq
     DO jk = 1 , klev , 1
      zqx ( jk , jm ) = pclv ( jl , jk , jm ) + ptsphy * tendency_tmp_cld ( jl , jk , jm )
      zqx0 ( jk , jm ) = pclv ( jl , jk , jm ) + ptsphy * tendency_tmp_cld ( jl , jk , jm )
     END DO
    END DO
!$acc loop seq
    DO jm = 1 , 5 , 1
!$acc loop seq
     DO jk = 1 , klev + 1 , 1
      zpfplsx ( jk , jm ) = 0.0_jprb
     END DO
    END DO
!$acc loop seq
    DO jm = 1 , 5 , 1
!$acc loop seq
     DO jk = 1 , klev , 1
      zqxn2d ( jk , jm ) = 0.0_jprb
      zlneg ( jk , jm ) = 0.0_jprb
     END DO
    END DO
    prainfrac_toprfz ( jl ) = 0.0_jprb
!$acc loop seq
    DO jk = 1 , klev , 1
     IF ( zqx ( jk , 1 ) + zqx ( jk , 2 ) < rlmin .OR. za ( jk ) < ramin ) THEN
      zlneg ( jk , 1 ) = zlneg ( jk , 1 ) + zqx ( jk , 1 )
      zqadj = zqx ( jk , 1 ) * zqtmst
      tendency_loc_q ( jl , jk ) = tendency_loc_q ( jl , jk ) + zqadj
      tendency_loc_t ( jl , jk ) = tendency_loc_t ( jl , jk ) - ralvdcp * zqadj
      zqx ( jk , 5 ) = zqx ( jk , 5 ) + zqx ( jk , 1 )
      zqx ( jk , 1 ) = 0.0_jprb
      zlneg ( jk , 2 ) = zlneg ( jk , 2 ) + zqx ( jk , 2 )
      zqadj = zqx ( jk , 2 ) * zqtmst
      tendency_loc_q ( jl , jk ) = tendency_loc_q ( jl , jk ) + zqadj
      tendency_loc_t ( jl , jk ) = tendency_loc_t ( jl , jk ) - ralsdcp * zqadj
      zqx ( jk , 5 ) = zqx ( jk , 5 ) + zqx ( jk , 2 )
      zqx ( jk , 2 ) = 0.0_jprb
      za ( jk ) = 0.0_jprb
     END IF
    END DO
!$acc loop seq
    DO jm = 1 , 4 , 1
!$acc loop seq
     DO jk = 1 , klev , 1
      IF ( zqx ( jk , jm ) < rlmin ) THEN
       zlneg ( jk , jm ) = zlneg ( jk , jm ) + zqx ( jk , jm )
       zqadj = zqx ( jk , jm ) * zqtmst
       tendency_loc_q ( jl , jk ) = tendency_loc_q ( jl , jk ) + zqadj
       IF ( iphase ( jm ) == 1 ) THEN
        tendency_loc_t ( jl , jk ) = tendency_loc_t ( jl , jk ) - ralvdcp * zqadj
       END IF
       IF ( iphase ( jm ) == 2 ) THEN
        tendency_loc_t ( jl , jk ) = tendency_loc_t ( jl , jk ) - ralsdcp * zqadj
       END IF
       zqx ( jk , 5 ) = zqx ( jk , 5 ) + zqx ( jk , jm )
       zqx ( jk , jm ) = 0.0_jprb
      END IF
     END DO
    END DO
!$acc loop seq
    DO jk = 1 , klev , 1
     zfoealfa ( jk ) = real ( min ( 1.0_jprb , ( ( max ( rtice , min ( rtwat , ztp1 ( jk ) ) ) - rtice ) * rtwat_rtice_r ) ** ( 2&
      ) ) , kind =JPRB )
     zfoeewmt ( jk ) = min ( real ( r2es * ( real ( min ( 1.0_jprb , ( ( max ( rtice , min ( rtwat , ztp1 ( jk ) ) ) - rtice ) *&
      rtwat_rtice_r ) ** ( 2 ) ) , kind =JPRB ) * exp ( r3les * ( ztp1 ( jk ) - rtt ) / ( ztp1 ( jk )&
      - r4les ) ) + ( 1.0_jprb - real ( min ( 1.0_jprb , ( ( max ( rtice , min ( rtwat , ztp1 ( jk ) ) ) - rtice ) * rtwat_rtice_r&
      ) ** ( 2 ) ) , kind =JPRB ) ) * exp ( r3ies * ( ztp1 ( jk ) - rtt ) / ( ztp1 ( jk ) - r4ies ) )&
      ) , kind =JPRB ) / pap ( jl , jk ) , 0.5_jprb )
     zqsmix ( jk ) = zfoeewmt ( jk )
     zqsmix ( jk ) = zqsmix ( jk ) / ( 1.0_jprb - retv * zqsmix ( jk ) )
     zalfa = real ( max ( 0.0_jprb , sign ( 1.0_jprb , ztp1 ( jk ) - rtt ) ) , kind =JPRB )
     zfoeew ( jk ) = min ( ( zalfa * real ( r2es * exp ( r3les * ( ztp1 ( jk ) - rtt ) / ( ztp1 ( jk ) - r4les ) ) , kind =&
     JPRB ) + ( 1.0_jprb - zalfa ) * real ( r2es * exp ( r3ies * ( ztp1 ( jk ) - rtt ) / ( ztp1 ( jk )&
      - r4ies ) ) , kind =JPRB ) ) / pap ( jl , jk ) , 0.5_jprb )
     zfoeew ( jk ) = min ( 0.5_jprb , zfoeew ( jk ) )
     zqsice ( jk ) = zfoeew ( jk ) / ( 1.0_jprb - retv * zfoeew ( jk ) )
     zfoeeliqt ( jk ) = min ( real ( r2es * exp ( r3les * ( ztp1 ( jk ) - rtt ) / ( ztp1 ( jk ) - r4les ) ) , kind =&
     JPRB ) / pap ( jl , jk ) , 0.5_jprb )
     zqsliq ( jk ) = zfoeeliqt ( jk )
     zqsliq ( jk ) = zqsliq ( jk ) / ( 1.0_jprb - retv * zqsliq ( jk ) )
    END DO
!$acc loop seq
    DO jk = 1 , klev , 1
     za ( jk ) = max ( 0.0_jprb , min ( 1.0_jprb , za ( jk ) ) )
     zli ( jk ) = zqx ( jk , 1 ) + zqx ( jk , 2 )
     IF ( zli ( jk ) > rlmin ) THEN
      zliqfrac ( jk ) = zqx ( jk , 1 ) / zli ( jk )
      zicefrac ( jk ) = 1.0_jprb - zliqfrac ( jk )
     ELSE
      zliqfrac ( jk ) = 0.0_jprb
      zicefrac ( jk ) = 0.0_jprb
     END IF
    END DO
    ztrpaus ( jl ) = 0.1_jprb
    zpaphd ( jl ) = 1.0_jprb / paph ( jl , klev + 1 )
!$acc loop seq
    DO jk = 1 , klev - 1 , 1
     zsig = pap ( jl , jk ) * zpaphd ( jl )
     IF ( zsig > 0.1_jprb .AND. ztp1 ( jk ) > ztp1 ( jk + 1 ) .AND. zsig < 0.4_jprb ) THEN
      ztrpaus ( jl ) = zsig
     END IF
    END DO
    zanewm1 ( jl ) = 0.0_jprb
    zda ( jl ) = 0.0_jprb
    zcovpclr ( jl ) = 0.0_jprb
    zcovpmax ( jl ) = 0.0_jprb
    zcovptot ( jl ) = 0.0_jprb
    zcldtopdist ( jl ) = 0.0_jprb
!$acc loop seq
    DO jk = ncldtop , klev , 1
!$acc loop seq
     DO jm = 1 , 5 , 1
      zqxfg ( jm ) = zqx ( jk , jm )
     END DO
     zlicld ( jl ) = 0.0_jprb
     zrainaut ( jl ) = 0.0_jprb
     zrainacc ( jl ) = 0.0_jprb
     zsnowaut ( jl ) = 0.0_jprb
     zldefr ( jl ) = 0.0_jprb
     zacust ( jl ) = 0.0_jprb
     zqpretot ( jl ) = 0.0_jprb
     zlfinalsum ( jl ) = 0.0_jprb
     zlcond1 ( jl ) = 0.0_jprb
     zlcond2 ( jl ) = 0.0_jprb
     zsupsat ( jl ) = 0.0_jprb
     zlevapl ( jl ) = 0.0_jprb
     zlevapi ( jl ) = 0.0_jprb
     zsolab ( jl ) = 0.0_jprb
     zsolac ( jl ) = 0.0_jprb
     zicetot ( jl ) = 0.0_jprb
!$acc loop seq
     DO jm = 1 , 5 , 1
!$acc loop seq
      DO jn = 1 , 5 , 1
       zsolqb ( jn , jm ) = 0.0_jprb
       zsolqa ( jn , jm ) = 0.0_jprb
      END DO
     END DO
!$acc loop seq
     DO jm = 1 , 5 , 1
      zfallsrce ( jm ) = 0.0_jprb
      zfallsink ( jm ) = 0.0_jprb
      zconvsrce ( jm ) = 0.0_jprb
      zconvsink ( jm ) = 0.0_jprb
      zpsupsatsrce ( jm ) = 0.0_jprb
      zratio ( jm ) = 0.0_jprb
     END DO
     IF ( llcldbudcc ) THEN
      zbudcc ( : ) = 0.0_jprb
     END IF
     IF ( llcldbudl ) THEN
      zbudl ( : ) = 0.0_jprb
     END IF
     IF ( llcldbudi ) THEN
      zbudi ( : ) = 0.0_jprb
     END IF
     zdp ( jl ) = paph ( jl , jk + 1 ) - paph ( jl , jk )
     zgdp ( jl ) = rg / zdp ( jl )
     zrho ( jl ) = pap ( jl , jk ) / ( rd * ztp1 ( jk ) )
     zdtgdp ( jl ) = ptsphy * zgdp ( jl )
     zrdtgdp ( jl ) = zdp ( jl ) * ( 1.0_jprb / ( ptsphy * rg ) )
     IF ( jk > 1 ) THEN
      zdtgdpf ( jl ) = ptsphy * rg / ( pap ( jl , jk ) - pap ( jl , jk - 1 ) )
     END IF
     zfacw = r5les / ( ztp1 ( jk ) - r4les ) ** ( 2 )
     zcor = 1.0_jprb / ( 1.0_jprb - retv * zfoeeliqt ( jk ) )
     zdqsliqdt ( jl ) = zfacw * zcor * zqsliq ( jk )
     zcorqsliq ( jl ) = 1.0_jprb + ralvdcp * zdqsliqdt ( jl )
     zfaci = r5ies / ( ztp1 ( jk ) - r4ies ) ** ( 2 )
     zcor = 1.0_jprb / ( 1.0_jprb - retv * zfoeew ( jk ) )
     zdqsicedt ( jl ) = zfaci * zcor * zqsice ( jk )
     zcorqsice ( jl ) = 1.0_jprb + ralsdcp * zdqsicedt ( jl )
     zalfaw = zfoealfa ( jk )
     zalfawm ( jl ) = zalfaw
     zfac = zalfaw * zfacw + ( 1.0_jprb - zalfaw ) * zfaci
     zcor = 1.0_jprb / ( 1.0_jprb - retv * zfoeewmt ( jk ) )
     zdqsmixdt ( jl ) = zfac * zcor * zqsmix ( jk )
     zcorqsmix ( jl ) = 1.0_jprb + real ( real ( min ( 1.0_jprb , ( ( max ( rtice , min ( rtwat , ztp1 ( jk ) ) ) - rtice ) *&
      rtwat_rtice_r ) ** ( 2 ) ) , kind =JPRB ) * ralvdcp + ( 1.0_jprb - real ( min ( 1.0_jprb , ( (&
      max ( rtice , min ( rtwat , ztp1 ( jk ) ) ) - rtice ) * rtwat_rtice_r ) ** ( 2 ) ) , kind = JPRB &
      ) ) * ralsdcp , kind =JPRB ) * zdqsmixdt ( jl )
     zevaplimmix ( jl ) = max ( ( zqsmix ( jk ) - zqx ( jk , 5 ) ) / zcorqsmix ( jl ) , 0.0_jprb )
     zevaplimliq ( jl ) = max ( ( zqsliq ( jk ) - zqx ( jk , 5 ) ) / zcorqsliq ( jl ) , 0.0_jprb )
     zevaplimice ( jl ) = max ( ( zqsice ( jk ) - zqx ( jk , 5 ) ) / zcorqsice ( jl ) , 0.0_jprb )
     ztmpa = 1.0_jprb * 1d0 / max ( za ( jk ) , zepsec )
     zliqcld ( jl ) = zqx ( jk , 1 ) * ztmpa
     zicecld ( jl ) = zqx ( jk , 2 ) * ztmpa
     zlicld ( jl ) = zliqcld ( jl ) + zicecld ( jl )
     IF ( zqx ( jk , 1 ) < rlmin ) THEN
      zsolqa ( 5 , 1 ) = zqx ( jk , 1 )
      zsolqa ( 1 , 5 ) = ( - zqx ( jk , 1 ) )
     END IF
     IF ( zqx ( jk , 2 ) < rlmin ) THEN
      zsolqa ( 5 , 2 ) = zqx ( jk , 2 )
      zsolqa ( 2 , 5 ) = ( - zqx ( jk , 2 ) )
     END IF
     zfokoop ( jl ) = real ( min ( rkoop1 - rkoop2 * ztp1 ( jk ) , real ( r2es * exp ( r3les * ( ztp1 ( jk ) - rtt ) / ( ztp1 ( jk&
      ) - r4les ) ) , kind =JPRB ) * 1_JPRB / real ( r2es * exp ( r3ies * ( ztp1 ( jk ) - rtt ) / ( ztp1&
      ( jk ) - r4ies ) ) , kind =JPRB ) ) , kind =JPRB )
     IF ( nssopt == 0 .OR. ztp1 ( jk ) >= rtt ) THEN
      zfac = 1.0_jprb
      zfaci = 1.0_jprb
     ELSE
      zfac = za ( jk ) + zfokoop ( jl ) * ( 1.0_jprb - za ( jk ) )
      zfaci = ptsphy / rkooptau
     END IF
     IF ( za ( jk ) > 1.0_jprb - ramin ) THEN
      zsupsat ( jl ) = max ( ( zqx ( jk , 5 ) - zfac * zqsice ( jk ) ) / zcorqsice ( jl ) , 0.0_jprb )
     ELSE
      zqp1env = ( zqx ( jk , 5 ) - za ( jk ) * zqsice ( jk ) ) * 1d0 / max ( 1.0_jprb - za ( jk ) , zepsilon )
      zsupsat ( jl ) = max ( ( 1.0_jprb - za ( jk ) ) * ( zqp1env - zfac * zqsice ( jk ) ) / zcorqsice ( jl ) , 0.0_jprb )
     END IF
     IF ( zsupsat ( jl ) > zepsec ) THEN
      IF ( ztp1 ( jk ) > rthomo ) THEN
       zsolqa ( 1 , 5 ) = zsolqa ( 1 , 5 ) + zsupsat ( jl )
       zsolqa ( 5 , 1 ) = zsolqa ( 5 , 1 ) - zsupsat ( jl )
       zqxfg ( 1 ) = zqxfg ( 1 ) + zsupsat ( jl )
      ELSE
       zsolqa ( 2 , 5 ) = zsolqa ( 2 , 5 ) + zsupsat ( jl )
       zsolqa ( 5 , 2 ) = zsolqa ( 5 , 2 ) - zsupsat ( jl )
       zqxfg ( 2 ) = zqxfg ( 2 ) + zsupsat ( jl )
      END IF
      zsolac ( jl ) = ( 1.0_jprb - za ( jk ) ) * zfaci
      IF ( llcldbudl .AND. ztp1 ( jk ) > rthomo ) THEN
       zbudl ( 1 ) = zsupsat ( jl ) * zqtmst
      END IF
      IF ( llcldbudi .AND. ztp1 ( jk ) <= rthomo ) THEN
       zbudi ( 1 ) = zsupsat ( jl ) * zqtmst
      END IF
      IF ( llcldbudcc ) THEN
       zbudcc ( 1 ) = zsolac ( jl ) * zqtmst
      END IF
     END IF
     IF ( psupsat ( jl , jk ) > zepsec ) THEN
      IF ( ztp1 ( jk ) > rthomo ) THEN
       zsolqa ( 1 , 1 ) = zsolqa ( 1 , 1 ) + psupsat ( jl , jk )
       zpsupsatsrce ( 1 ) = psupsat ( jl , jk )
       zqxfg ( 1 ) = zqxfg ( 1 ) + psupsat ( jl , jk )
       IF ( llcldbudl ) THEN
        zbudl ( 2 ) = psupsat ( jl , jk ) * zqtmst
       END IF
      ELSE
       zsolqa ( 2 , 2 ) = zsolqa ( 2 , 2 ) + psupsat ( jl , jk )
       zpsupsatsrce ( 2 ) = psupsat ( jl , jk )
       zqxfg ( 2 ) = zqxfg ( 2 ) + psupsat ( jl , jk )
       IF ( llcldbudi ) THEN
        zbudi ( 2 ) = psupsat ( jl , jk ) * zqtmst
       END IF
      END IF
      zsolac ( jl ) = ( 1.0_jprb - za ( jk ) ) * zfaci
      IF ( llcldbudcc ) THEN
       zbudcc ( 2 ) = zsolac ( jl ) * zqtmst
      END IF
     END IF
     IF ( jk >= ncldtop .AND. jk < klev ) THEN
      plude ( jl , jk ) = plude ( jl , jk ) * zdtgdp ( jl )
      IF ( ldcum ( jl ) .AND. plu ( jl , jk + 1 ) > zepsec .AND. plude ( jl , jk ) > rlmin ) THEN
       zsolac ( jl ) = zsolac ( jl ) + plude ( jl , jk ) / plu ( jl , jk + 1 )
       zalfaw = zfoealfa ( jk )
       zconvsrce ( 1 ) = zalfaw * plude ( jl , jk )
       zconvsrce ( 2 ) = ( 1.0_jprb - zalfaw ) * plude ( jl , jk )
       zsolqa ( 1 , 1 ) = zsolqa ( 1 , 1 ) + zconvsrce ( 1 )
       zsolqa ( 2 , 2 ) = zsolqa ( 2 , 2 ) + zconvsrce ( 2 )
       IF ( llcldbudl ) THEN
        zbudl ( 3 ) = zconvsrce ( 1 ) * zqtmst
       END IF
       IF ( llcldbudi ) THEN
        zbudi ( 3 ) = zconvsrce ( 2 ) * zqtmst
       END IF
       IF ( llcldbudcc ) THEN
        zbudcc ( 3 ) = zqtmst * plude ( jl , jk ) / plu ( jl , jk + 1 )
       END IF
      ELSE
       plude ( jl , jk ) = 0.0_jprb
      END IF
      IF ( ldcum ( jl ) ) THEN
       zsolqa ( 4 , 4 ) = zsolqa ( 4 , 4 ) + psnde ( jl , jk ) * zdtgdp ( jl )
      END IF
     END IF
     IF ( jk > ncldtop ) THEN
      zmf ( jl ) = max ( 0.0_jprb , ( pmfu ( jl , jk ) + pmfd ( jl , jk ) ) * zdtgdp ( jl ) )
      zacust ( jl ) = zmf ( jl ) * zanewm1 ( jl )
!$acc loop seq
      DO jm = 1 , 5 , 1
       IF ( ( .NOT. llfall ( jm ) ) .AND. iphase ( jm ) > 0 ) THEN
        zlcust ( jm ) = zmf ( jl ) * zqxnm1 ( jm )
        zconvsrce ( jm ) = zconvsrce ( jm ) + zlcust ( jm )
       END IF
      END DO
      zdtdp = zrdcp * 0.5_jprb * ( ztp1 ( jk - 1 ) + ztp1 ( jk ) ) / paph ( jl , jk )
      zdtforc = zdtdp * ( pap ( jl , jk ) - pap ( jl , jk - 1 ) )
      zdqs ( jl ) = zanewm1 ( jl ) * zdtforc * zdqsmixdt ( jl )
!$acc loop seq
      DO jm = 1 , 5 , 1
       IF ( ( .NOT. llfall ( jm ) ) .AND. iphase ( jm ) > 0 ) THEN
        zlfinal = max ( 0.0_jprb , zlcust ( jm ) - zdqs ( jl ) )
        zevap = min ( zlcust ( jm ) - zlfinal , zevaplimmix ( jl ) )
        zlfinal = zlcust ( jm ) - zevap
        zlfinalsum ( jl ) = zlfinalsum ( jl ) + zlfinal
        zsolqa ( jm , jm ) = zsolqa ( jm , jm ) + zlcust ( jm )
        zsolqa ( 5 , jm ) = zsolqa ( 5 , jm ) + zevap
        zsolqa ( jm , 5 ) = zsolqa ( jm , 5 ) - zevap
        IF ( llcldbudl .AND. jm == 1 ) THEN
        zbudl ( 4 ) = zlcust ( jm ) * zqtmst
        END IF
        IF ( llcldbudi .AND. jm == 2 ) THEN
        zbudi ( 4 ) = zlcust ( jm ) * zqtmst
        END IF
        IF ( llcldbudl .AND. jm == 1 ) THEN
        zbudl ( 5 ) = ( - zevap * zqtmst )
        END IF
        IF ( llcldbudi .AND. jm == 2 ) THEN
        zbudi ( 5 ) = ( - zevap * zqtmst )
        END IF
       END IF
      END DO
      IF ( zlfinalsum ( jl ) < zepsec ) THEN
       zacust ( jl ) = 0.0_jprb
      END IF
      zsolac ( jl ) = zsolac ( jl ) + zacust ( jl )
      IF ( llcldbudcc ) THEN
       zbudcc ( 4 ) = zacust ( jl ) * zqtmst
      END IF
     END IF
     IF ( jk < klev ) THEN
      zmfdn = max ( 0.0_jprb , ( pmfu ( jl , jk + 1 ) + pmfd ( jl , jk + 1 ) ) * zdtgdp ( jl ) )
      zsolab ( jl ) = zsolab ( jl ) + zmfdn
      zsolqb ( 1 , 1 ) = zsolqb ( 1 , 1 ) + zmfdn
      zsolqb ( 2 , 2 ) = zsolqb ( 2 , 2 ) + zmfdn
      zconvsink ( 1 ) = zmfdn
      zconvsink ( 2 ) = zmfdn
     END IF
     zldifdt ( jl ) = rcldiff * ptsphy
     IF ( ktype ( jl ) > 0 .AND. plude ( jl , jk ) > zepsec ) THEN
      zldifdt ( jl ) = rcldiff_convi * zldifdt ( jl )
     END IF
     IF ( zli ( jk ) > zepsec ) THEN
      ze = zldifdt ( jl ) * max ( zqsmix ( jk ) - zqx ( jk , 5 ) , 0.0_jprb )
      zleros = za ( jk ) * ze
      zleros = min ( zleros , zevaplimmix ( jl ) )
      zleros = min ( zleros , zli ( jk ) )
      zaeros = zleros / zlicld ( jl )
      zsolac ( jl ) = zsolac ( jl ) - zaeros
      zsolqa ( 5 , 1 ) = zsolqa ( 5 , 1 ) + zliqfrac ( jk ) * zleros
      zsolqa ( 1 , 5 ) = zsolqa ( 1 , 5 ) - zliqfrac ( jk ) * zleros
      zsolqa ( 5 , 2 ) = zsolqa ( 5 , 2 ) + zicefrac ( jk ) * zleros
      zsolqa ( 2 , 5 ) = zsolqa ( 2 , 5 ) - zicefrac ( jk ) * zleros
      IF ( llcldbudl ) THEN
       zbudl ( 7 ) = ( - zliqfrac ( jk ) * zleros * zqtmst )
      END IF
      IF ( llcldbudi ) THEN
       zbudi ( 7 ) = ( - zicefrac ( jk ) * zleros * zqtmst )
      END IF
      IF ( llcldbudcc ) THEN
       zbudcc ( 7 ) = ( - zaeros * zqtmst )
      END IF
     END IF
     zdtdp = zrdcp * ztp1 ( jk ) / pap ( jl , jk )
     zdpmxdt = zdp ( jl ) * zqtmst
     zmfdn = 0.0_jprb
     IF ( jk < klev ) THEN
      zmfdn = pmfu ( jl , jk + 1 ) + pmfd ( jl , jk + 1 )
     END IF
     zwtot = pvervel ( jl , jk ) + 0.5_jprb * rg * ( pmfu ( jl , jk ) + pmfd ( jl , jk ) + zmfdn )
     zwtot = min ( zdpmxdt , max ( ( - zdpmxdt ) , zwtot ) )
     zzzdt = phrsw ( jl , jk ) + phrlw ( jl , jk )
     zdtdiab = min ( zdpmxdt * zdtdp , max ( ( - zdpmxdt * zdtdp ) , zzzdt ) ) * ptsphy + ralfdcp * zldefr ( jl )
     zdtforc = zdtdp * zwtot * ptsphy + zdtdiab
     zqold ( jl ) = zqsmix ( jk )
     ztold ( jl ) = ztp1 ( jk )
     ztp1 ( jk ) = ztp1 ( jk ) + zdtforc
     ztp1 ( jk ) = max ( ztp1 ( jk ) , 160.0_jprb )
     llflag = .TRUE.
     zqp = 1.0_jprb / pap ( jl , jk )
     zqsat = real ( r2es * ( real ( min ( 1.0_jprb , ( ( max ( rtice , min ( rtwat , ztp1 ( jk ) ) ) - rtice ) * rtwat_rtice_r )&
      ** ( 2 ) ) , kind =JPRB ) * exp ( r3les * ( ztp1 ( jk ) - rtt ) / ( ztp1 ( jk ) - r4les ) ) + (&
      1.0_jprb - real ( min ( 1.0_jprb , ( ( max ( rtice , min ( rtwat , ztp1 ( jk ) ) ) - rtice ) * rtwat_rtice_r ) ** ( 2 ) ) ,&
      kind =JPRB ) ) * exp ( r3ies * ( ztp1 ( jk ) - rtt ) / ( ztp1 ( jk ) - r4ies ) ) ) , kind =&
     JPRB ) * zqp
     zqsat = min ( 0.5_jprb , zqsat )
     zcor = 1.0_jprb / ( 1.0_jprb - retv * zqsat )
     zqsat = zqsat * zcor
     zcond = ( zqsmix ( jk ) - zqsat ) * 1d0 / ( 1.0_jprb + zqsat * zcor * real ( real ( min ( 1.0_jprb , ( ( max ( rtice , min (&
      rtwat , ztp1 ( jk ) ) ) - rtice ) * rtwat_rtice_r ) ** ( 2 ) ) , kind =JPRB ) * r5alvcp * (&
      1.0_jprb / ( ztp1 ( jk ) - r4les ) ** ( 2 ) ) + ( 1.0_jprb - real ( min ( 1.0_jprb , ( ( max ( rtice , min ( rtwat , ztp1 (&
      jk ) ) ) - rtice ) * rtwat_rtice_r ) ** ( 2 ) ) , kind =JPRB ) ) * r5alscp * ( 1.0_jprb / ( ztp1&
      ( jk ) - r4ies ) ** ( 2 ) ) , kind =JPRB ) )
     ztp1 ( jk ) = ztp1 ( jk ) + real ( real ( min ( 1.0_jprb , ( ( max ( rtice , min ( rtwat , ztp1 ( jk ) ) ) - rtice ) *&
      rtwat_rtice_r ) ** ( 2 ) ) , kind =JPRB ) * ralvdcp + ( 1.0_jprb - real ( min ( 1.0_jprb , ( (&
      max ( rtice , min ( rtwat , ztp1 ( jk ) ) ) - rtice ) * rtwat_rtice_r ) ** ( 2 ) ) , kind = JPRB &
      ) ) * ralsdcp , kind =JPRB ) * zcond
     zqsmix ( jk ) = zqsmix ( jk ) - zcond
     zqsat = real ( r2es * ( real ( min ( 1.0_jprb , ( ( max ( rtice , min ( rtwat , ztp1 ( jk ) ) ) - rtice ) * rtwat_rtice_r )&
      ** ( 2 ) ) , kind =JPRB ) * exp ( r3les * ( ztp1 ( jk ) - rtt ) / ( ztp1 ( jk ) - r4les ) ) + (&
      1.0_jprb - real ( min ( 1.0_jprb , ( ( max ( rtice , min ( rtwat , ztp1 ( jk ) ) ) - rtice ) * rtwat_rtice_r ) ** ( 2 ) ) ,&
      kind =JPRB ) ) * exp ( r3ies * ( ztp1 ( jk ) - rtt ) / ( ztp1 ( jk ) - r4ies ) ) ) , kind =&
     JPRB ) * zqp
     zqsat = min ( 0.5_jprb , zqsat )
     zcor = 1.0_jprb / ( 1.0_jprb - retv * zqsat )
     zqsat = zqsat * zcor
     zcond1 = ( zqsmix ( jk ) - zqsat ) * 1d0 / ( 1.0_jprb + zqsat * zcor * real ( real ( min ( 1.0_jprb , ( ( max ( rtice , min (&
      rtwat , ztp1 ( jk ) ) ) - rtice ) * rtwat_rtice_r ) ** ( 2 ) ) , kind =JPRB ) * r5alvcp * (&
      1.0_jprb / ( ztp1 ( jk ) - r4les ) ** ( 2 ) ) + ( 1.0_jprb - real ( min ( 1.0_jprb , ( ( max ( rtice , min ( rtwat , ztp1 (&
      jk ) ) ) - rtice ) * rtwat_rtice_r ) ** ( 2 ) ) , kind =JPRB ) ) * r5alscp * ( 1.0_jprb / ( ztp1&
      ( jk ) - r4ies ) ** ( 2 ) ) , kind =JPRB ) )
     ztp1 ( jk ) = ztp1 ( jk ) + real ( real ( min ( 1.0_jprb , ( ( max ( rtice , min ( rtwat , ztp1 ( jk ) ) ) - rtice ) *&
      rtwat_rtice_r ) ** ( 2 ) ) , kind =JPRB ) * ralvdcp + ( 1.0_jprb - real ( min ( 1.0_jprb , ( (&
      max ( rtice , min ( rtwat , ztp1 ( jk ) ) ) - rtice ) * rtwat_rtice_r ) ** ( 2 ) ) , kind = JPRB &
      ) ) * ralsdcp , kind =JPRB ) * zcond1
     zqsmix ( jk ) = zqsmix ( jk ) - zcond1
     zdqs ( jl ) = zqsmix ( jk ) - zqold ( jl )
     zqsmix ( jk ) = zqold ( jl )
     ztp1 ( jk ) = ztold ( jl )
     IF ( zdqs ( jl ) > 0.0_jprb ) THEN
      zlevap = za ( jk ) * min ( zdqs ( jl ) , zlicld ( jl ) )
      zlevap = min ( zlevap , zevaplimmix ( jl ) )
      zlevap = min ( zlevap , max ( zqsmix ( jk ) - zqx ( jk , 5 ) , 0.0_jprb ) )
      zlevapl ( jl ) = zliqfrac ( jk ) * zlevap
      zlevapi ( jl ) = zicefrac ( jk ) * zlevap
      zsolqa ( 5 , 1 ) = zsolqa ( 5 , 1 ) + zliqfrac ( jk ) * zlevap
      zsolqa ( 1 , 5 ) = zsolqa ( 1 , 5 ) - zliqfrac ( jk ) * zlevap
      zsolqa ( 5 , 2 ) = zsolqa ( 5 , 2 ) + zicefrac ( jk ) * zlevap
      zsolqa ( 2 , 5 ) = zsolqa ( 2 , 5 ) - zicefrac ( jk ) * zlevap
      IF ( llcldbudl ) THEN
       zbudl ( 8 ) = ( - zliqfrac ( jk ) * zlevap * zqtmst )
      END IF
      IF ( llcldbudi ) THEN
       zbudi ( 8 ) = ( - zicefrac ( jk ) * zlevap * zqtmst )
      END IF
     END IF
     IF ( zdqs ( jl ) <= ( - rlmin ) .AND. za ( jk ) > zepsec ) THEN
      zlcond1 ( jl ) = max ( ( - zdqs ( jl ) ) , 0.0_jprb )
      IF ( za ( jk ) > 0.99_jprb ) THEN
       zcor = 1.0_jprb / ( 1.0_jprb - retv * zqsmix ( jk ) )
       zcdmax = ( zqx ( jk , 5 ) - zqsmix ( jk ) ) * 1d0 / ( 1.0_jprb + zcor * zqsmix ( jk ) * real ( real ( min ( 1.0_jprb , ( (&
        max ( rtice , min ( rtwat , ztp1 ( jk ) ) ) - rtice ) * rtwat_rtice_r ) ** ( 2 ) ) , kind = JPRB &
        ) * r5alvcp * ( 1.0_jprb / ( ztp1 ( jk ) - r4les ) ** ( 2 ) ) + ( 1.0_jprb - real ( min ( 1.0_jprb , ( ( max ( rtice ,&
        min ( rtwat , ztp1 ( jk ) ) ) - rtice ) * rtwat_rtice_r ) ** ( 2 ) ) , kind =JPRB ) ) *&
        r5alscp * ( 1.0_jprb / ( ztp1 ( jk ) - r4ies ) ** ( 2 ) ) , kind =JPRB ) )
      ELSE
       zcdmax = ( zqx ( jk , 5 ) - za ( jk ) * zqsmix ( jk ) ) / za ( jk )
      END IF
      zlcond1 ( jl ) = max ( min ( zlcond1 ( jl ) , zcdmax ) , 0.0_jprb )
      zlcond1 ( jl ) = za ( jk ) * zlcond1 ( jl )
      IF ( zlcond1 ( jl ) < rlmin ) THEN
       zlcond1 ( jl ) = 0.0_jprb
      END IF
      IF ( ztp1 ( jk ) > rthomo ) THEN
       zsolqa ( 1 , 5 ) = zsolqa ( 1 , 5 ) + zlcond1 ( jl )
       zsolqa ( 5 , 1 ) = zsolqa ( 5 , 1 ) - zlcond1 ( jl )
       zqxfg ( 1 ) = zqxfg ( 1 ) + zlcond1 ( jl )
       IF ( llcldbudl ) THEN
        zbudl ( 9 ) = zlcond1 ( jl ) * zqtmst
       END IF
      ELSE
       zsolqa ( 2 , 5 ) = zsolqa ( 2 , 5 ) + zlcond1 ( jl )
       zsolqa ( 5 , 2 ) = zsolqa ( 5 , 2 ) - zlcond1 ( jl )
       zqxfg ( 2 ) = zqxfg ( 2 ) + zlcond1 ( jl )
       IF ( llcldbudi ) THEN
        zbudi ( 9 ) = zlcond1 ( jl ) * zqtmst
       END IF
      END IF
     END IF
     IF ( zdqs ( jl ) <= ( - rlmin ) .AND. za ( jk ) < 1.0_jprb - zepsec ) THEN
      zrhc = ramid
      zsigk = pap ( jl , jk ) / paph ( jl , klev + 1 )
      IF ( zsigk > 0.8_jprb ) THEN
       zrhc = ramid + ( 1.0_jprb - ramid ) * ( ( zsigk - 0.8_jprb ) / 0.2_jprb ) ** ( 2 )
      END IF
      IF ( nssopt == 0 ) THEN
       zqe = ( zqx ( jk , 5 ) - za ( jk ) * zqsice ( jk ) ) * 1d0 / max ( zepsec , 1.0_jprb - za ( jk ) )
       zqe = max ( 0.0_jprb , zqe )
      ELSE
       IF ( nssopt == 1 ) THEN
        zqe = ( zqx ( jk , 5 ) - za ( jk ) * zqsice ( jk ) ) * 1d0 / max ( zepsec , 1.0_jprb - za ( jk ) )
        zqe = max ( 0.0_jprb , zqe )
       ELSE
        IF ( nssopt == 2 ) THEN
        zqe = zqx ( jk , 5 )
        ELSE
        IF ( nssopt == 3 ) THEN
        zqe = zqx ( jk , 5 ) + zli ( jk )
        END IF
        END IF
       END IF
      END IF
      IF ( nssopt == 0 .OR. ztp1 ( jk ) >= rtt ) THEN
       zfac = 1.0_jprb
      ELSE
       zfac = zfokoop ( jl )
      END IF
      IF ( zqe >= zrhc * zqsice ( jk ) * zfac .AND. zqe < zqsice ( jk ) * zfac ) THEN
       zacond = ( - ( 1.0_jprb - za ( jk ) ) * zfac * zdqs ( jl ) * 1d0 / max ( 2.0_jprb * ( zfac * zqsice ( jk ) - zqe ) , zepsec&
        ) )
       zacond = min ( zacond , 1.0_jprb - za ( jk ) )
       zlcond2 ( jl ) = ( - zfac * zdqs ( jl ) * 0.5_jprb * zacond )
       zzdl = 2.0_jprb * ( zfac * zqsice ( jk ) - zqe ) * 1d0 / max ( zepsec , 1.0_jprb - za ( jk ) )
       IF ( zfac * zdqs ( jl ) < ( - zzdl ) ) THEN
        zlcondlim = ( za ( jk ) - 1.0_jprb ) * zfac * zdqs ( jl ) - zfac * zqsice ( jk ) + zqx ( jk , 5 )
        zlcond2 ( jl ) = min ( zlcond2 ( jl ) , zlcondlim )
       END IF
       zlcond2 ( jl ) = max ( zlcond2 ( jl ) , 0.0_jprb )
       IF ( zlcond2 ( jl ) < rlmin .OR. 1.0_jprb - za ( jk ) < zepsec ) THEN
        zlcond2 ( jl ) = 0.0_jprb
        zacond = 0.0_jprb
       END IF
       IF ( zlcond2 ( jl ) == 0.0_jprb ) THEN
        zacond = 0.0_jprb
       END IF
       zsolac ( jl ) = zsolac ( jl ) + zacond
       IF ( llcldbudcc ) THEN
        zbudcc ( 10 ) = zacond * zqtmst
       END IF
       IF ( ztp1 ( jk ) > rthomo ) THEN
        zsolqa ( 1 , 5 ) = zsolqa ( 1 , 5 ) + zlcond2 ( jl )
        zsolqa ( 5 , 1 ) = zsolqa ( 5 , 1 ) - zlcond2 ( jl )
        zqxfg ( 1 ) = zqxfg ( 1 ) + zlcond2 ( jl )
        IF ( llcldbudl ) THEN
        zbudl ( 10 ) = zlcond2 ( jl ) * zqtmst
        END IF
       ELSE
        zsolqa ( 2 , 5 ) = zsolqa ( 2 , 5 ) + zlcond2 ( jl )
        zsolqa ( 5 , 2 ) = zsolqa ( 5 , 2 ) - zlcond2 ( jl )
        zqxfg ( 2 ) = zqxfg ( 2 ) + zlcond2 ( jl )
        IF ( llcldbudi ) THEN
        zbudi ( 10 ) = zlcond2 ( jl ) * zqtmst
        END IF
       END IF
      END IF
     END IF
     IF ( idepice == 1 ) THEN
      IF ( za ( jk ) >= rcldtopcf .AND. za ( jk - 1 ) < rcldtopcf ) THEN
       zcldtopdist ( jl ) = 0.0_jprb
      ELSE
       zcldtopdist ( jl ) = zcldtopdist ( jl ) + zdp ( jl ) / ( zrho ( jl ) * rg )
      END IF
      IF ( zqxfg ( 1 ) > rlmin .AND. ztp1 ( jk ) < rtt ) THEN
       zvpice = real ( r2es * exp ( r3ies * ( ztp1 ( jk ) - rtt ) / ( ztp1 ( jk ) - r4ies ) ) , kind = JPRB &
        ) * rv / rd
       zvpliq = zvpice * zfokoop ( jl )
       zicenuclei ( jl ) = 1000.0_jprb * exp ( 12.96_jprb * ( zvpliq - zvpice ) / zvpliq - 0.639_jprb )
       zadd = rlstt * ( rlstt / ( rv * ztp1 ( jk ) ) - 1.0_jprb ) / ( 0.024_jprb * ztp1 ( jk ) )
       zbdd = rv * ztp1 ( jk ) * pap ( jl , jk ) / ( 2.21_jprb * zvpice )
       zcvds = 7.8_jprb * ( zicenuclei ( jl ) / zrho ( jl ) ) ** ( 0.666_jprb ) * ( zvpliq - zvpice ) / ( 8.87_jprb * ( zadd +&
        zbdd ) * zvpice )
       zice0 = max ( zicecld ( jl ) , zicenuclei ( jl ) * riceinit / zrho ( jl ) )
       zinew = ( 0.666_jprb * zcvds * ptsphy + zice0 ** ( 0.666_jprb ) ) ** ( 1.5_jprb )
       zdepos = max ( za ( jk ) * ( zinew - zice0 ) , 0.0_jprb )
       zdepos = min ( zdepos , zqxfg ( 1 ) )
       zinfactor = min ( zicenuclei ( jl ) / 15000.0_jprb , 1.0_jprb )
       zdepos = zdepos * min ( zinfactor + ( 1.0_jprb - zinfactor ) * ( rdepliqrefrate + zcldtopdist ( jl ) / rdepliqrefdepth ) ,&
        1.0_jprb )
       zsolqa ( 2 , 1 ) = zsolqa ( 2 , 1 ) + zdepos
       zsolqa ( 1 , 2 ) = zsolqa ( 1 , 2 ) - zdepos
       zqxfg ( 2 ) = zqxfg ( 2 ) + zdepos
       zqxfg ( 1 ) = zqxfg ( 1 ) - zdepos
       IF ( llcldbudl ) THEN
        zbudl ( 11 ) = ( - zdepos * zqtmst )
       END IF
       IF ( llcldbudi ) THEN
        zbudi ( 11 ) = zdepos * zqtmst
       END IF
      END IF
     ELSE
      IF ( idepice == 2 ) THEN
       IF ( za ( jk ) >= rcldtopcf .AND. za ( jk - 1 ) < rcldtopcf ) THEN
        zcldtopdist ( jl ) = 0.0_jprb
       ELSE
        zcldtopdist ( jl ) = zcldtopdist ( jl ) + zdp ( jl ) / ( zrho ( jl ) * rg )
       END IF
       IF ( zqxfg ( 1 ) > rlmin .AND. ztp1 ( jk ) < rtt ) THEN
        zvpice = real ( r2es * exp ( r3ies * ( ztp1 ( jk ) - rtt ) / ( ztp1 ( jk ) - r4ies ) ) , kind = JPRB &
         ) * rv / rd
        zvpliq = zvpice * zfokoop ( jl )
        zicenuclei ( jl ) = 1000.0_jprb * exp ( 12.96_jprb * ( zvpliq - zvpice ) / zvpliq - 0.639_jprb )
        zice0 = max ( zicecld ( jl ) , zicenuclei ( jl ) * riceinit / zrho ( jl ) )
        ztcg = 1.0_jprb
        zfacx1i = 1.0_jprb
        zaplusb = rcl_apb1 * zvpice - rcl_apb2 * zvpice * ztp1 ( jk ) + pap ( jl , jk ) * rcl_apb3 * ztp1 ( jk ) ** ( 3.0_jprb )
        zcorrfac = sqrt ( 1.0_jprb / zrho ( jl ) )
        zcorrfac2 = ( ztp1 ( jk ) / 273.0_jprb ) ** ( 1.5_jprb ) * ( 393.0_jprb / ( ztp1 ( jk ) + 120.0_jprb ) )
        zpr02 = zrho ( jl ) * zice0 * rcl_const1i / ( ztcg * zfacx1i )
        zterm1 = ( zvpliq - zvpice ) * ztp1 ( jk ) ** ( 2.0_jprb ) * zvpice * zcorrfac2 * ztcg * rcl_const2i * zfacx1i / ( zrho (&
         jl ) * zaplusb * zvpice )
        zterm2 = 0.65_jprb * rcl_const6i * zpr02 ** ( rcl_const4i ) + rcl_const3i * sqrt ( zcorrfac ) * sqrt ( zrho ( jl ) ) *&
         zpr02 ** ( rcl_const5i ) / sqrt ( zcorrfac2 )
        zdepos = max ( za ( jk ) * zterm1 * zterm2 * ptsphy , 0.0_jprb )
        zdepos = min ( zdepos , zqxfg ( 1 ) )
        zinfactor = min ( zicenuclei ( jl ) / 15000.0_jprb , 1.0_jprb )
        zdepos = zdepos * min ( zinfactor + ( 1.0_jprb - zinfactor ) * ( rdepliqrefrate + zcldtopdist ( jl ) / rdepliqrefdepth ) ,&
         1.0_jprb )
        zsolqa ( 2 , 1 ) = zsolqa ( 2 , 1 ) + zdepos
        zsolqa ( 1 , 2 ) = zsolqa ( 1 , 2 ) - zdepos
        zqxfg ( 2 ) = zqxfg ( 2 ) + zdepos
        zqxfg ( 1 ) = zqxfg ( 1 ) - zdepos
        IF ( llcldbudl ) THEN
        zbudl ( 11 ) = ( - zdepos * zqtmst )
        END IF
        IF ( llcldbudi ) THEN
        zbudi ( 11 ) = zdepos * zqtmst
        END IF
       END IF
      END IF
     END IF
     ztmpa = 1.0_jprb * 1d0 / max ( za ( jk ) , zepsec )
     zliqcld ( jl ) = zqxfg ( 1 ) * ztmpa
     zicecld ( jl ) = zqxfg ( 2 ) * ztmpa
     zlicld ( jl ) = zliqcld ( jl ) + zicecld ( jl )
!$acc loop seq
     DO jm = 1 , 5 , 1
      IF ( llfall ( jm ) .OR. jm == 2 ) THEN
       IF ( jk > ncldtop ) THEN
        zfallsrce ( jm ) = zpfplsx ( jk , jm ) * zdtgdp ( jl )
        zsolqa ( jm , jm ) = zsolqa ( jm , jm ) + zfallsrce ( jm )
        zqxfg ( jm ) = zqxfg ( jm ) + zfallsrce ( jm )
        zqpretot ( jl ) = zqpretot ( jl ) + zqxfg ( jm )
        IF ( llcldbudi .AND. jm == 2 ) THEN
        zbudi ( 12 ) = zfallsrce ( jm ) * zqtmst
        END IF
       END IF
       IF ( laericesed .AND. jm == 2 ) THEN
        zre_ice = pre_ice ( jl , jk )
        zvqx ( 2 ) = 0.002_jprb * zre_ice ** ( 1.0_jprb )
       END IF
       zfall = zvqx ( jm ) * zrho ( jl )
       zfallsink ( jm ) = zdtgdp ( jl ) * zfall
      END IF
     END DO
     IF ( zqpretot ( jl ) > zepsec ) THEN
      zcovptot ( jl ) = 1.0_jprb - ( 1.0_jprb - zcovptot ( jl ) ) * ( 1.0_jprb - max ( za ( jk ) , za ( jk - 1 ) ) ) * 1d0 / (&
       1.0_jprb - min ( za ( jk - 1 ) , 1.0_jprb - 1.0e-6_jprb ) )
      zcovptot ( jl ) = max ( zcovptot ( jl ) , rcovpmin )
      zcovpclr ( jl ) = max ( 0.0_jprb , zcovptot ( jl ) - za ( jk ) )
      zraincld ( jl ) = zqxfg ( 3 ) / zcovptot ( jl )
      zsnowcld ( jl ) = zqxfg ( 4 ) / zcovptot ( jl )
      zcovpmax ( jl ) = max ( zcovptot ( jl ) , zcovpmax ( jl ) )
     ELSE
      zraincld ( jl ) = 0.0_jprb
      zsnowcld ( jl ) = 0.0_jprb
      zcovptot ( jl ) = 0.0_jprb
      zcovpclr ( jl ) = 0.0_jprb
      zcovpmax ( jl ) = 0.0_jprb
     END IF
     IF ( ztp1 ( jk ) <= rtt ) THEN
      IF ( zicecld ( jl ) > zepsec ) THEN
       zzco = ptsphy * rsnowlin1 * exp ( rsnowlin2 * ( ztp1 ( jk ) - rtt ) )
       IF ( laericeauto ) THEN
        zlcrit = picrit_aer ( jl , jk )
        zzco = zzco * ( rnice / pnice ( jl , jk ) ) ** ( 0.333_jprb )
       ELSE
        zlcrit = rlcritsnow
       END IF
       zsnowaut ( jl ) = zzco * ( 1.0_jprb - exp ( ( - ( zicecld ( jl ) / zlcrit ) ** ( 2 ) ) ) )
       zsolqb ( 4 , 2 ) = zsolqb ( 4 , 2 ) + zsnowaut ( jl )
      END IF
     END IF
     IF ( zliqcld ( jl ) > zepsec ) THEN
      IF ( iwarmrain == 1 ) THEN
       zzco = rkconv * ptsphy
       IF ( laerliqautolsp ) THEN
        zlcrit = plcrit_aer ( jl , jk )
        zzco = zzco * ( rccn / pccn ( jl , jk ) ) ** ( 0.333_jprb )
       ELSE
        IF ( plsm ( jl ) > 0.5_jprb ) THEN
        zlcrit = rclcrit_land
        ELSE
        zlcrit = rclcrit_sea
        END IF
       END IF
       zprecip = ( zpfplsx ( jk , 4 ) + zpfplsx ( jk , 3 ) ) * 1d0 / max ( zepsec , zcovptot ( jl ) )
       zcfpr = 1.0_jprb + rprc1 * sqrt ( max ( zprecip , 0.0_jprb ) )
       IF ( laerliqcoll ) THEN
        zcfpr = zcfpr * ( rccn / pccn ( jl , jk ) ) ** ( 0.333_jprb )
       END IF
       zzco = zzco * zcfpr
       zlcrit = zlcrit * 1d0 / max ( zcfpr , zepsec )
       IF ( zliqcld ( jl ) / zlcrit < 20.0_jprb ) THEN
        zrainaut ( jl ) = zzco * ( 1.0_jprb - exp ( ( - ( zliqcld ( jl ) / zlcrit ) ** ( 2 ) ) ) )
       ELSE
        zrainaut ( jl ) = zzco
       END IF
       IF ( ztp1 ( jk ) <= rtt ) THEN
        zsolqb ( 4 , 1 ) = zsolqb ( 4 , 1 ) + zrainaut ( jl )
       ELSE
        zsolqb ( 3 , 1 ) = zsolqb ( 3 , 1 ) + zrainaut ( jl )
       END IF
      ELSE
       IF ( iwarmrain == 2 ) THEN
        IF ( plsm ( jl ) > 0.5_jprb ) THEN
        zconst = rcl_kk_cloud_num_land
        zlcrit = rclcrit_land
        ELSE
        zconst = rcl_kk_cloud_num_sea
        zlcrit = rclcrit_sea
        END IF
        IF ( zliqcld ( jl ) > zlcrit ) THEN
        zrainaut ( jl ) = 1.5_jprb * za ( jk ) * ptsphy * rcl_kkaau * zliqcld ( jl ) ** ( rcl_kkbauq ) * zconst ** ( rcl_kkbaun )
        zrainaut ( jl ) = min ( zrainaut ( jl ) , zqxfg ( 1 ) )
        IF ( zrainaut ( jl ) < zepsec ) THEN
        zrainaut ( jl ) = 0.0_jprb
        END IF
        zrainacc ( jl ) = 2.0_jprb * za ( jk ) * ptsphy * rcl_kkaac * ( zliqcld ( jl ) * zraincld ( jl ) ) ** ( rcl_kkbac )
        zrainacc ( jl ) = min ( zrainacc ( jl ) , zqxfg ( 1 ) )
        IF ( zrainacc ( jl ) < zepsec ) THEN
        zrainacc ( jl ) = 0.0_jprb
        END IF
        ELSE
        zrainaut ( jl ) = 0.0_jprb
        zrainacc ( jl ) = 0.0_jprb
        END IF
        IF ( ztp1 ( jk ) <= rtt ) THEN
        zsolqa ( 4 , 1 ) = zsolqa ( 4 , 1 ) + zrainaut ( jl )
        zsolqa ( 4 , 1 ) = zsolqa ( 4 , 1 ) + zrainacc ( jl )
        zsolqa ( 1 , 4 ) = zsolqa ( 1 , 4 ) - zrainaut ( jl )
        zsolqa ( 1 , 4 ) = zsolqa ( 1 , 4 ) - zrainacc ( jl )
        ELSE
        zsolqa ( 3 , 1 ) = zsolqa ( 3 , 1 ) + zrainaut ( jl )
        zsolqa ( 3 , 1 ) = zsolqa ( 3 , 1 ) + zrainacc ( jl )
        zsolqa ( 1 , 3 ) = zsolqa ( 1 , 3 ) - zrainaut ( jl )
        zsolqa ( 1 , 3 ) = zsolqa ( 1 , 3 ) - zrainacc ( jl )
        END IF
       END IF
      END IF
     END IF
     IF ( iwarmrain > 1 ) THEN
      IF ( ztp1 ( jk ) <= rtt .AND. zliqcld ( jl ) > zepsec ) THEN
       zfallcorr = ( rdensref / zrho ( jl ) ) ** ( 0.4_jprb )
       IF ( zcovptot ( jl ) > 0.01_jprb .AND. zsnowcld ( jl ) > zepsec ) THEN
        zsnowrime ( jl ) = 0.3_jprb * zcovptot ( jl ) * ptsphy * rcl_const7s * zfallcorr * ( zrho ( jl ) * zsnowcld ( jl ) *&
         rcl_const1s ) ** ( rcl_const8s )
        zsnowrime ( jl ) = min ( zsnowrime ( jl ) , 1.0_jprb )
        zsolqb ( 4 , 1 ) = zsolqb ( 4 , 1 ) + zsnowrime ( jl )
       END IF
      END IF
     END IF
     zicetot ( jl ) = zqxfg ( 2 ) + zqxfg ( 4 )
     zmeltmax ( jl ) = 0.0_jprb
     IF ( zicetot ( jl ) > zepsec .AND. ztp1 ( jk ) > rtt ) THEN
      zsubsat = max ( zqsice ( jk ) - zqx ( jk , 5 ) , 0.0_jprb )
      ztdmtw0 = ztp1 ( jk ) - rtt - zsubsat * ( ztw1 + ztw2 * ( pap ( jl , jk ) - ztw3 ) - ztw4 * ( ztp1 ( jk ) - ztw5 ) )
      zcons1 = abs ( ptsphy * ( 1.0_jprb + 0.5_jprb * ztdmtw0 ) / rtaumel )
      zmeltmax ( jl ) = max ( ztdmtw0 * zcons1 * zrldcp , 0.0_jprb )
     END IF
!$acc loop seq
     DO jm = 1 , 5 , 1
      IF ( iphase ( jm ) == 2 ) THEN
       jn = imelt ( jm )
       IF ( zicetot ( jl ) > zepsec .AND. zmeltmax ( jl ) > zepsec ) THEN
        zalfa = zqxfg ( jm ) / zicetot ( jl )
        zmelt = min ( zqxfg ( jm ) , zalfa * zmeltmax ( jl ) )
        zqxfg ( jm ) = zqxfg ( jm ) - zmelt
        zqxfg ( jn ) = zqxfg ( jn ) + zmelt
        zsolqa ( jn , jm ) = zsolqa ( jn , jm ) + zmelt
        zsolqa ( jm , jn ) = zsolqa ( jm , jn ) - zmelt
        IF ( llcldbudi .AND. jm == 2 ) THEN
        zbudi ( 15 ) = ( - zmelt * zqtmst )
        END IF
        IF ( llcldbudi .AND. jm == 4 ) THEN
        zbudi ( 16 ) = ( - zmelt * zqtmst )
        END IF
        IF ( llcldbudl .AND. jm == 2 ) THEN
        zbudl ( 17 ) = zmelt * zqtmst
        END IF
       END IF
      END IF
     END DO
     IF ( zqx ( jk , 3 ) > zepsec ) THEN
      IF ( ztp1 ( jk ) <= rtt .AND. ztp1 ( jk - 1 ) > rtt ) THEN
       zqpretot ( jl ) = max ( zqx ( jk , 4 ) + zqx ( jk , 3 ) , zepsec )
       prainfrac_toprfz ( jl ) = zqx ( jk , 3 ) / zqpretot ( jl )
      END IF
      IF ( ztp1 ( jk ) < rtt ) THEN
       IF ( prainfrac_toprfz ( jl ) > 0.8 ) THEN
        zlambda = ( rcl_fac1 / ( zrho ( jl ) * zqx ( jk , 3 ) ) ) ** ( rcl_fac2 )
        ztemp = rcl_fzrab * ( ztp1 ( jk ) - rtt )
        zfrz = ptsphy * ( rcl_const5r / zrho ( jl ) ) * ( exp ( ztemp ) - 1.0_jprb ) * zlambda ** ( rcl_const6r )
        zfrzmax ( jl ) = max ( zfrz , 0.0_jprb )
       ELSE
        zcons1 = abs ( ptsphy * ( 1.0_jprb + 0.5_jprb * ( rtt - ztp1 ( jk ) ) ) / rtaumel )
        zfrzmax ( jl ) = max ( ( rtt - ztp1 ( jk ) ) * zcons1 * zrldcp , 0.0_jprb )
       END IF
       IF ( zfrzmax ( jl ) > zepsec ) THEN
        zfrz = min ( zqx ( jk , 3 ) , zfrzmax ( jl ) )
        zsolqa ( 4 , 3 ) = zsolqa ( 4 , 3 ) + zfrz
        zsolqa ( 3 , 4 ) = zsolqa ( 3 , 4 ) - zfrz
        IF ( llcldbudl ) THEN
        zbudl ( 18 ) = zfrz * zqtmst
        END IF
       END IF
      END IF
     END IF
     zfrzmax ( jl ) = max ( ( rthomo - ztp1 ( jk ) ) * zrldcp , 0.0_jprb )
     jm = 1
     jn = imelt ( jm )
     IF ( zfrzmax ( jl ) > zepsec .AND. zqxfg ( jm ) > zepsec ) THEN
      zfrz = min ( zqxfg ( jm ) , zfrzmax ( jl ) )
      zsolqa ( jn , jm ) = zsolqa ( jn , jm ) + zfrz
      zsolqa ( jm , jn ) = zsolqa ( jm , jn ) - zfrz
      IF ( llcldbudl ) THEN
       zbudl ( 19 ) = zfrz * zqtmst
      END IF
     END IF
     IF ( ievaprain == 1 ) THEN
      zzrh = rprecrhmax + ( 1.0_jprb - rprecrhmax ) * zcovpmax ( jl ) * 1d0 / max ( zepsec , 1.0_jprb - za ( jk ) )
      zzrh = min ( max ( zzrh , rprecrhmax ) , 1.0_jprb )
      zqe = ( zqx ( jk , 5 ) - za ( jk ) * zqsliq ( jk ) ) * 1d0 / max ( zepsec , 1.0_jprb - za ( jk ) )
      zqe = max ( 0.0_jprb , min ( zqe , zqsliq ( jk ) ) )
      llo1 = zcovpclr ( jl ) > zepsec .AND. zqxfg ( 3 ) > zepsec .AND. zqe < zzrh * zqsliq ( jk )
      IF ( llo1 ) THEN
       zpreclr = zqxfg ( 3 ) * zcovpclr ( jl ) * 1d0 / sign ( max ( abs ( zcovptot ( jl ) * zdtgdp ( jl ) ) , zepsilon ) ,&
        zcovptot ( jl ) * zdtgdp ( jl ) )
       zbeta1 = sqrt ( pap ( jl , jk ) / paph ( jl , klev + 1 ) ) / rvrfactor * zpreclr * 1d0 / max ( zcovpclr ( jl ) , zepsec )
       zbeta = rg * rpecons * 0.5_jprb * zbeta1 ** ( 0.5777_jprb )
       zdenom = 1.0_jprb + zbeta * ptsphy * zcorqsliq ( jl )
       zdpr = zcovpclr ( jl ) * zbeta * ( zqsliq ( jk ) - zqe ) / zdenom * zdp ( jl ) * zrg_r
       zdpevap = zdpr * zdtgdp ( jl )
       zevap = min ( zdpevap , zqxfg ( 3 ) )
       zsolqa ( 5 , 3 ) = zsolqa ( 5 , 3 ) + zevap
       zsolqa ( 3 , 5 ) = zsolqa ( 3 , 5 ) - zevap
       zcovptot ( jl ) = max ( rcovpmin , zcovptot ( jl ) - max ( 0.0_jprb , ( zcovptot ( jl ) - za ( jk ) ) * zevap / zqxfg ( 3 )&
        ) )
       zqxfg ( 3 ) = zqxfg ( 3 ) - zevap
      END IF
     ELSE
      IF ( ievaprain == 2 ) THEN
       zzrh = rprecrhmax + ( 1.0_jprb - rprecrhmax ) * zcovpmax ( jl ) * 1d0 / max ( zepsec , 1.0_jprb - za ( jk ) )
       zzrh = min ( max ( zzrh , rprecrhmax ) , 1.0_jprb )
       zzrh = min ( 0.8_jprb , zzrh )
       zqe = max ( 0.0_jprb , min ( zqx ( jk , 5 ) , zqsliq ( jk ) ) )
       llo1 = zcovpclr ( jl ) > zepsec .AND. zqxfg ( 3 ) > zepsec .AND. zqe < zzrh * zqsliq ( jk )
       IF ( llo1 ) THEN
        zpreclr = zqxfg ( 3 ) / zcovptot ( jl )
        zfallcorr = ( rdensref / zrho ( jl ) ) ** ( 0.4_jprb )
        zesatliq = rv / rd * real ( r2es * exp ( r3les * ( ztp1 ( jk ) - rtt ) / ( ztp1 ( jk ) - r4les ) ) , kind =&
        JPRB )
        zlambda = ( rcl_fac1 / ( zrho ( jl ) * zpreclr ) ) ** ( rcl_fac2 )
        zevap_denom = rcl_cdenom1 * zesatliq - rcl_cdenom2 * ztp1 ( jk ) * zesatliq + rcl_cdenom3 * ztp1 ( jk ) ** ( 3.0_jprb ) *&
         pap ( jl , jk )
        zcorr2 = ( ztp1 ( jk ) / 273.0_jprb ) ** ( 1.5_jprb ) * 393.0_jprb / ( ztp1 ( jk ) + 120.0_jprb )
        zka = rcl_ka273 * zcorr2
        zsubsat = max ( zzrh * zqsliq ( jk ) - zqe , 0.0_jprb )
        zbeta = 0.5_jprb / zqsliq ( jk ) * ztp1 ( jk ) ** ( 2.0_jprb ) * zesatliq * rcl_const1r * ( zcorr2 / zevap_denom ) * (&
         0.78_jprb / zlambda ** ( rcl_const4r ) + rcl_const2r * sqrt ( zrho ( jl ) * zfallcorr ) / ( sqrt ( zcorr2 ) * zlambda **&
         ( rcl_const3r ) ) )
        zdenom = 1.0_jprb + zbeta * ptsphy
        zdpevap = zcovpclr ( jl ) * zbeta * ptsphy * zsubsat / zdenom
        zevap = min ( zdpevap , zqxfg ( 3 ) )
        zsolqa ( 5 , 3 ) = zsolqa ( 5 , 3 ) + zevap
        zsolqa ( 3 , 5 ) = zsolqa ( 3 , 5 ) - zevap
        zcovptot ( jl ) = max ( rcovpmin , zcovptot ( jl ) - max ( 0.0_jprb , ( zcovptot ( jl ) - za ( jk ) ) * zevap / zqxfg ( 3&
         ) ) )
        zqxfg ( 3 ) = zqxfg ( 3 ) - zevap
       END IF
      END IF
     END IF
     IF ( ievapsnow == 1 ) THEN
      zzrh = rprecrhmax + ( 1.0_jprb - rprecrhmax ) * zcovpmax ( jl ) * 1d0 / max ( zepsec , 1.0_jprb - za ( jk ) )
      zzrh = min ( max ( zzrh , rprecrhmax ) , 1.0_jprb )
      zqe = ( zqx ( jk , 5 ) - za ( jk ) * zqsice ( jk ) ) * 1d0 / max ( zepsec , 1.0_jprb - za ( jk ) )
      zqe = max ( 0.0_jprb , min ( zqe , zqsice ( jk ) ) )
      llo1 = zcovpclr ( jl ) > zepsec .AND. zqxfg ( 4 ) > zepsec .AND. zqe < zzrh * zqsice ( jk )
      IF ( llo1 ) THEN
       zpreclr = zqxfg ( 4 ) * zcovpclr ( jl ) * 1d0 / sign ( max ( abs ( zcovptot ( jl ) * zdtgdp ( jl ) ) , zepsilon ) ,&
        zcovptot ( jl ) * zdtgdp ( jl ) )
       zbeta1 = sqrt ( pap ( jl , jk ) / paph ( jl , klev + 1 ) ) / rvrfactor * zpreclr * 1d0 / max ( zcovpclr ( jl ) , zepsec )
       zbeta = rg * rpecons * zbeta1 ** ( 0.5777_jprb )
       zdenom = 1.0_jprb + zbeta * ptsphy * zcorqsice ( jl )
       zdpr = zcovpclr ( jl ) * zbeta * ( zqsice ( jk ) - zqe ) / zdenom * zdp ( jl ) * zrg_r
       zdpevap = zdpr * zdtgdp ( jl )
       zevap = min ( zdpevap , zqxfg ( 4 ) )
       zsolqa ( 5 , 4 ) = zsolqa ( 5 , 4 ) + zevap
       zsolqa ( 4 , 5 ) = zsolqa ( 4 , 5 ) - zevap
       zcovptot ( jl ) = max ( rcovpmin , zcovptot ( jl ) - max ( 0.0_jprb , ( zcovptot ( jl ) - za ( jk ) ) * zevap / zqxfg ( 4 )&
        ) )
       zqxfg ( 4 ) = zqxfg ( 4 ) - zevap
      END IF
     ELSE
      IF ( ievapsnow == 2 ) THEN
       zzrh = rprecrhmax + ( 1.0_jprb - rprecrhmax ) * zcovpmax ( jl ) * 1d0 / max ( zepsec , 1.0_jprb - za ( jk ) )
       zzrh = min ( max ( zzrh , rprecrhmax ) , 1.0_jprb )
       zqe = ( zqx ( jk , 5 ) - za ( jk ) * zqsice ( jk ) ) * 1d0 / max ( zepsec , 1.0_jprb - za ( jk ) )
       zqe = max ( 0.0_jprb , min ( zqe , zqsice ( jk ) ) )
       llo1 = zcovpclr ( jl ) > zepsec .AND. zqx ( jk , 4 ) > zepsec .AND. zqe < zzrh * zqsice ( jk )
       IF ( llo1 ) THEN
        zpreclr = zqx ( jk , 4 ) / zcovptot ( jl )
        zvpice = real ( r2es * exp ( r3ies * ( ztp1 ( jk ) - rtt ) / ( ztp1 ( jk ) - r4ies ) ) , kind = JPRB &
         ) * rv / rd
        ztcg = 1.0_jprb
        zfacx1s = 1.0_jprb
        zaplusb = rcl_apb1 * zvpice - rcl_apb2 * zvpice * ztp1 ( jk ) + pap ( jl , jk ) * rcl_apb3 * ztp1 ( jk ) ** ( 3 )
        zcorrfac = sqrt ( 1.0_jprb / zrho ( jl ) )
        zcorrfac2 = ( ztp1 ( jk ) / 273.0_jprb ) ** ( 1.5_jprb ) * ( 393.0_jprb / ( ztp1 ( jk ) + 120.0_jprb ) )
        zpr02 = zrho ( jl ) * zpreclr * rcl_const1s / ( ztcg * zfacx1s )
        zterm1 = ( zqsice ( jk ) - zqe ) * ztp1 ( jk ) ** ( 2 ) * zvpice * zcorrfac2 * ztcg * rcl_const2s * zfacx1s / ( zrho ( jl&
         ) * zaplusb * zqsice ( jk ) )
        zterm2 = 0.65 * rcl_const6s * zpr02 ** ( rcl_const4s ) + rcl_const3s * sqrt ( zcorrfac ) * sqrt ( zrho ( jl ) ) * zpr02 **&
         ( rcl_const5s ) / sqrt ( zcorrfac2 )
        zdpevap = max ( zcovpclr ( jl ) * zterm1 * zterm2 * ptsphy , 0.0_jprb )
        zevap = min ( zdpevap , zevaplimice ( jl ) )
        zevap = min ( zevap , zqx ( jk , 4 ) )
        zsolqa ( 5 , 4 ) = zsolqa ( 5 , 4 ) + zevap
        zsolqa ( 4 , 5 ) = zsolqa ( 4 , 5 ) - zevap
        zcovptot ( jl ) = max ( rcovpmin , zcovptot ( jl ) - max ( 0.0_jprb , ( zcovptot ( jl ) - za ( jk ) ) * zevap / zqx ( jk ,&
         4 ) ) )
        zqxfg ( 4 ) = zqxfg ( 4 ) - zevap
       END IF
      END IF
     END IF
!$acc loop seq
     DO jm = 1 , 5 , 1
      IF ( llfall ( jm ) ) THEN
       IF ( zqxfg ( jm ) < rlmin ) THEN
        zsolqa ( 5 , jm ) = zsolqa ( 5 , jm ) + zqxfg ( jm )
        zsolqa ( jm , 5 ) = zsolqa ( jm , 5 ) - zqxfg ( jm )
       END IF
      END IF
     END DO
     zanew = ( za ( jk ) + zsolac ( jl ) ) / ( 1.0_jprb + zsolab ( jl ) )
     zanew = min ( zanew , 1.0_jprb )
     IF ( zanew < ramin ) THEN
      zanew = 0.0_jprb
     END IF
     zda ( jl ) = zanew - zaorig ( jk )
     zanewm1 ( jl ) = zanew
!$acc loop seq
     DO jm = 1 , 5 , 1
      DO jn = 1 , 5 , 1
       llindex3 ( jn , jm ) = .FALSE.
      END DO
      zsinksum ( jm ) = 0.0_jprb
     END DO
!$acc loop seq
     DO jm = 1 , 5 , 1
!$acc loop seq
      DO jn = 1 , 5 , 1
       zsinksum ( jm ) = zsinksum ( jm ) - zsolqa ( jm , jn )
      END DO
     END DO
!$acc loop seq
     DO jm = 1 , 5 , 1
      zmax = max ( zqx ( jk , jm ) , zepsec )
      zrat = max ( zsinksum ( jm ) , zmax )
      zratio ( jm ) = zmax / zrat
     END DO
!$acc loop seq
     DO jm = 1 , 5 , 1
      zsinksum ( jm ) = 0.0_jprb
     END DO
!$acc loop seq
     DO jm = 1 , 5 , 1
      psum_solqa ( jl ) = 0.0
!$acc loop seq
      DO jn = 1 , 5 , 1
       psum_solqa ( jl ) = psum_solqa ( jl ) + zsolqa ( jm , jn )
      END DO
      zsinksum ( jm ) = zsinksum ( jm ) - psum_solqa ( jl )
      zmm = max ( zqx ( jk , jm ) , zepsec )
      zrr = max ( zsinksum ( jm ) , zmm )
      zratio ( jm ) = zmm / zrr
      zzratio = zratio ( jm )
!$acc loop seq
      DO jn = 1 , 5 , 1
       IF ( zsolqa ( jm , jn ) < 0.0_jprb ) THEN
        zsolqa ( jm , jn ) = zsolqa ( jm , jn ) * zzratio
        zsolqa ( jn , jm ) = zsolqa ( jn , jm ) * zzratio
       END IF
      END DO
     END DO
!$acc loop seq
     DO jm = 1 , 5 , 1
!$acc loop seq
      DO jn = 1 , 5 , 1
       IF ( jn == jm ) THEN
        zqlhs ( jn , jm ) = 1.0_jprb + zfallsink ( jm )
!$acc loop seq
        DO jo = 1 , 5 , 1
        zqlhs ( jn , jm ) = zqlhs ( jn , jm ) + zsolqb ( jo , jn )
        END DO
       ELSE
        zqlhs ( jn , jm ) = ( - zsolqb ( jn , jm ) )
       END IF
      END DO
     END DO
!$acc loop seq
     DO jm = 1 , 5 , 1
      zexplicit = 0.0_jprb
!$acc loop seq
      DO jn = 1 , 5 , 1
       zexplicit = zexplicit + zsolqa ( jm , jn )
      END DO
      zqxn ( jm ) = zqx ( jk , jm ) + zexplicit
     END DO
!$acc loop seq
     DO jn = 1 , 4 , 1
!$acc loop seq
      DO jm = jn + 1 , 5 , 1
       zqlhs ( jm , jn ) = zqlhs ( jm , jn ) / zqlhs ( jn , jn )
!$acc loop seq
       DO ik = jn + 1 , 5 , 1
        zqlhs ( jm , ik ) = zqlhs ( jm , ik ) - zqlhs ( jm , jn ) * zqlhs ( jn , ik )
       END DO
      END DO
     END DO
!$acc loop seq
     DO jn = 2 , 5 , 1
!$acc loop seq
      DO jm = 1 , jn - 1 , 1
       zqxn ( jn ) = zqxn ( jn ) - zqlhs ( jn , jm ) * zqxn ( jm )
      END DO
     END DO
     zqxn ( 5 ) = zqxn ( 5 ) / zqlhs ( 5 , 5 )
!$acc loop seq
     DO jn = 4 , 1 , (-1)
!$acc loop seq
      DO jm = jn + 1 , 5 , 1
       zqxn ( jn ) = zqxn ( jn ) - zqlhs ( jn , jm ) * zqxn ( jm )
      END DO
      zqxn ( jn ) = zqxn ( jn ) / zqlhs ( jn , jn )
     END DO
!$acc loop seq
     DO jn = 1 , 4 , 1
      IF ( zqxn ( jn ) < zepsec ) THEN
       zqxn ( 5 ) = zqxn ( 5 ) + zqxn ( jn )
       zqxn ( jn ) = 0.0_jprb
      END IF
     END DO
!$acc loop seq
     DO jm = 1 , 5 , 1
      zqxnm1 ( jm ) = zqxn ( jm )
      zqxn2d ( jk , jm ) = zqxn ( jm )
     END DO
!$acc loop seq
     DO jm = 1 , 5 , 1
      zpfplsx ( jk + 1 , jm ) = zfallsink ( jm ) * zqxn ( jm ) * zrdtgdp ( jl )
     END DO
     zqpretot ( jl ) = zpfplsx ( jk + 1 , 4 ) + zpfplsx ( jk + 1 , 3 )
     IF ( zqpretot ( jl ) < zepsec ) THEN
      zcovptot ( jl ) = 0.0_jprb
     END IF
!$acc loop seq
     DO jm = 1 , 4 , 1
      zfluxq ( jm ) = zpsupsatsrce ( jm ) + zconvsrce ( jm ) + zfallsrce ( jm ) - ( zfallsink ( jm ) + zconvsink ( jm ) ) * zqxn (&
       jm )
      IF ( iphase ( jm ) == 1 ) THEN
       tendency_loc_t ( jl , jk ) = tendency_loc_t ( jl , jk ) + ralvdcp * ( zqxn ( jm ) - zqx ( jk , jm ) - zfluxq ( jm ) ) *&
        zqtmst
      END IF
      IF ( iphase ( jm ) == 2 ) THEN
       tendency_loc_t ( jl , jk ) = tendency_loc_t ( jl , jk ) + ralsdcp * ( zqxn ( jm ) - zqx ( jk , jm ) - zfluxq ( jm ) ) *&
        zqtmst
      END IF
      tendency_loc_cld ( jl , jk , jm ) = tendency_loc_cld ( jl , jk , jm ) + ( zqxn ( jm ) - zqx0 ( jk , jm ) ) * zqtmst
     END DO
     tendency_loc_q ( jl , jk ) = tendency_loc_q ( jl , jk ) + ( zqxn ( 5 ) - zqx ( jk , 5 ) ) * zqtmst
     tendency_loc_a ( jl , jk ) = tendency_loc_a ( jl , jk ) + zda ( jl ) * zqtmst
     pcovptot ( jl , jk ) = zcovptot ( jl )
    END DO
!$acc loop seq
    DO jk = 1 , klev + 1 , 1
     pfplsl ( jl , jk ) = zpfplsx ( jk , 3 ) + zpfplsx ( jk , 1 )
     pfplsn ( jl , jk ) = zpfplsx ( jk , 4 ) + zpfplsx ( jk , 2 )
    END DO
    pfsqlf ( jl , 1 ) = 0.0_jprb
    pfsqif ( jl , 1 ) = 0.0_jprb
    pfsqrf ( jl , 1 ) = 0.0_jprb
    pfsqsf ( jl , 1 ) = 0.0_jprb
    pfcqlng ( jl , 1 ) = 0.0_jprb
    pfcqnng ( jl , 1 ) = 0.0_jprb
    pfcqrng ( jl , 1 ) = 0.0_jprb
    pfcqsng ( jl , 1 ) = 0.0_jprb
    pfsqltur ( jl , 1 ) = 0.0_jprb
    pfsqitur ( jl , 1 ) = 0.0_jprb
!$acc loop seq
    DO jk = 1 , klev , 1
     zgdph_r = ( - zrg_r * ( paph ( jl , jk + 1 ) - paph ( jl , jk ) ) * zqtmst )
     pfsqlf ( jl , jk + 1 ) = pfsqlf ( jl , jk )
     pfsqif ( jl , jk + 1 ) = pfsqif ( jl , jk )
     pfsqrf ( jl , jk + 1 ) = pfsqlf ( jl , jk )
     pfsqsf ( jl , jk + 1 ) = pfsqif ( jl , jk )
     pfcqlng ( jl , jk + 1 ) = pfcqlng ( jl , jk )
     pfcqnng ( jl , jk + 1 ) = pfcqnng ( jl , jk )
     pfcqrng ( jl , jk + 1 ) = pfcqlng ( jl , jk )
     pfcqsng ( jl , jk + 1 ) = pfcqnng ( jl , jk )
     pfsqltur ( jl , jk + 1 ) = pfsqltur ( jl , jk )
     pfsqitur ( jl , jk + 1 ) = pfsqitur ( jl , jk )
     zalfaw = zfoealfa ( jk )
     pfsqlf ( jl , jk + 1 ) = pfsqlf ( jl , jk + 1 ) + ( zqxn2d ( jk , 1 ) - zqx0 ( jk , 1 ) + pvfl ( jl , jk ) * ptsphy - zalfaw&
      * plude ( jl , jk ) ) * zgdph_r
     pfcqlng ( jl , jk + 1 ) = pfcqlng ( jl , jk + 1 ) + zlneg ( jk , 1 ) * zgdph_r
     pfsqltur ( jl , jk + 1 ) = pfsqltur ( jl , jk + 1 ) + pvfl ( jl , jk ) * ptsphy * zgdph_r
     pfsqrf ( jl , jk + 1 ) = pfsqrf ( jl , jk + 1 ) + ( zqxn2d ( jk , 3 ) - zqx0 ( jk , 3 ) ) * zgdph_r
     pfcqrng ( jl , jk + 1 ) = pfcqrng ( jl , jk + 1 ) + zlneg ( jk , 3 ) * zgdph_r
     pfsqif ( jl , jk + 1 ) = pfsqif ( jl , jk + 1 ) + ( zqxn2d ( jk , 2 ) - zqx0 ( jk , 2 ) + pvfi ( jl , jk ) * ptsphy - (&
      1.0_jprb - zalfaw ) * plude ( jl , jk ) ) * zgdph_r
     pfcqnng ( jl , jk + 1 ) = pfcqnng ( jl , jk + 1 ) + zlneg ( jk , 2 ) * zgdph_r
     pfsqitur ( jl , jk + 1 ) = pfsqitur ( jl , jk + 1 ) + pvfi ( jl , jk ) * ptsphy * zgdph_r
     pfsqsf ( jl , jk + 1 ) = pfsqsf ( jl , jk + 1 ) + ( zqxn2d ( jk , 4 ) - zqx0 ( jk , 4 ) ) * zgdph_r
     pfcqsng ( jl , jk + 1 ) = pfcqsng ( jl , jk + 1 ) + zlneg ( jk , 4 ) * zgdph_r
    END DO
!$acc loop seq
    DO jk = 1 , klev + 1 , 1
     pfhpsl ( jl , jk ) = ( - rlvtt * pfplsl ( jl , jk ) )
     pfhpsn ( jl , jk ) = ( - rlstt * pfplsn ( jl , jk ) )
    END DO
   END ASSOCIATE
  END DO

  end do

!$acc end parallel
!$acc end data

 END SUBROUTINE cloudsc_claw

END MODULE cloudsc_claw_mod

