! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

MODULE YOMCT0

USE PARKIND1  ,ONLY : JPIM, JPRB

IMPLICIT NONE

SAVE

!     ------------------------------------------------------------------

!*    Control variables for the job - constant within job
!     Values are identical for all models run under the OOPS layer.

!========== TECHNICAL SWITCHES ================================================

! ----- configuration:
! NCONF      : configuration of the job

!                0- 99 : 3-D integration job
!              100-199 : variational job
!              200-299 : 2-D integration job
!              300-349 : KALMAN filter
!              350-399 : predictability model             (currently unused)
!              400-499 : test of the adjoint
!              500-599 : test of the tangent linear model
!              600-699 : eigenvalue/vector solvers
!              700-799 : optimal interpolation
!              800-899 : sensitivity
!              900-999 : miscellaneous other configurations.

!                    1 : 3-D primitive equation model

!                  131 : incremental 4-D VAR/3-D VAR

!                  201 : shallow-water model
!                  202 : vorticity equation model

!                  302 : simplified extended Kalman filter (SEKF)

!                  401 : test of adjoint with 3-D P.E. model
!                  421 : test of adjoint with shallow-water model
!                  422 : test of adjoint with vorticity equation model

!                  501 : test of tangent linear with 3-D P.E. model
!                  521 : test of tangent linear with shallow-water model
!                  522 : test of tangent linear with vorticity equation model

!                  601 : eigenvalue/vector solver for 3-D P.E. model

!                  701 : optimal interpolation with CANARI

!                  801 : sensitivity with 3-D P.E. model

!                  901 : set up initial conditions (CPREP1)
!                  923 : initialisation of climatologic files
!                  931 : creation of an ARPEGE file containing the NESDIS SST (INCLITC).
!                  932 : interpolates the sea-ice concentration field from
!                        satellite data to the ARPEGE grid (CSEAICE).

! LR3D       : 3D model (currently NCONF in 1,131,302,401,501,601,801)
! LR2D       : 2D shallow water or vorticity model
! LRSHW      : 2D shallow water model
! LRVEQ      : 2D vorticity model

! ----- variables linked with diagnostics and outputs:
! LALLOPR    : .T. = print information about all allocations/deallocations

! ----- type of file used:
! LFDBOP     : .T. = fields data base utilized
! LGRBOP     : .T. = output in GRIB (not ARPEGE)
! LFBDAP     : .T. = diagnostics written on trajectory files
! LARPEGEF   : .T. = use ARPEGE files
! LARPEGEF_TRAJHR   : .T. = use ARPEGE files for high resolution trajectory
! LARPEGEF_TRAJBG   : .T. = use ARPEGE files for background
! LARPEGEF_RDGP_INIT     : .T. = use grid-point ARPEGE files
! LARPEGEF_RDGP_TRAJHR   : .T. = use grid-point ARPEGE files for HR trajectory
! LARPEGEF_RDGP_TRAJBG   : .T. = use grid-point ARPEGE files for background
! CNDISPP    : directory of display files

! ----- variables linked to post-processing:
! NFPOS      : = 0   <=> Fullpos inactive
! NFPOS      : = 1   <=> Fullpos active, standard backend post-processing
! NFPOS      : = 2   <=> Fullpos active, configuration for model geometry changes (new style)
! NFPOS      : = 927 <=> Fullpos active, configuration for model geometry changes (old style)
! LECFPOS    : if .TRUE., restrict the usage of Fullpos outputs to an alternative to model ones
! CFPNCF     : name of Full-POS control file (pseudo-configuration 927)
! LFPART2    : .T. = second part of interpolations for changing geometry
!              (pseudo-configuration 927)
! CFDIRLST   : path of postprocessing listing file
! CNPPATH    : directory of postprocessing namelist files

! ----- use of transmission coefficients stored in Fourier or spectral space:
! LRETCFOU   : .T. = reading Fourier transmission coefficients on file
!              for simplified physics.
! LWRTCFOU   : .T. = writing Fourier transmission coefficients on file.
! LRFOUTCNORM: activates diagnostics on Fourier transmission coefficients.
!              - LRFOUTCNORM(1)=.T.: global statistics.
!              - LRFOUTCNORM(2)=.T.: statistics per latitude.
!              - LRFOUTCNORM(3)=.T.: statistics per layer.
!              - LRFOUTCNORM(4)=.T.: statistics per wavenumber.
! LRGPTCNORM : activates diagnostics on grid-point transmission coefficients.
!              - LRGPTCNORM(1)=.T.: global statistics.
!              - LRGPTCNORM(2)=.T.: statistics per latitude.
!              - LRGPTCNORM(3)=.T.: statistics per layer.

! ----- other variables:
! LNF        : .T. = start, .F. = restart
! LSMSSIG    : .T. = send signals to SMS (ECMWF supervisor)
! CMETER     :  SMS or ECFLOW meter command (ECMWF supervisor)
! CEVENT     :  SMS or ECFLOW event command (ECMWF supervisor)
! LOPDIS     : .T. = calls OPDIS
! NCYCLE     : cycle identifier
! CNMEXP     : name of the experiment
!              An experiment is identified by its name (16 characters)
!              and its cycle (typically same experiment but without a bug)
! CFCLASS    : class for use by FDB
! CTYPE      : type for use by FDB
! LBACKG     : ???
! LSPBSBAL   : Read Jb Balance constraint files and apply Omega and non-linear
!              balance to SPBS vorticity perturbations
! LMINIM     : ???
! LIOLEVG    : .F. to enable the number of vertical levels in the model to be
!               less than the number of vertical levels in the files headers.
! NUNDEFLD   : index value for unused/undefined fields (default=-9999)
!            : should be set to 1 when using compiler subscript checking
! LREFOUT    : .T. compare to reference run
! LREFGEN    : .T. to generate reference file

!========== MODEL SWITCHES ====================================================
! Remark: this section containes some dynamical variables which cannot be
!         put in NAMDYN because they are used in routines called before SUDYN
!         and sometimes to allocate some arrays.

! ----- advection scheme, semi-Lagrangian scheme:
! LSLAG     : .TRUE. = semi-lagrangian on
! LTWOTL    : .TRUE. if two-time-level semi-Lagrangian scheme will be used

! ----- vertical discretisation, vertical boundaries:
! LREGETA   : .T.: for the interlayer L, ETA(L)=L/NFLEVG
!             .F.: for the interlayer L, ETA(L)=A(L)/P0+B(L)
! LVFE_REGETA: cf. LREGETA for "eta" used in VFE operators.

! ----- quadrature:
! NQUAD     : 1 ====> GAUSS

! ----- type of equations:
! LNHDYN    : .T. if 3-D non-hydrostatic dynamics is active

! ----- split ECMWF physics:
! LSLPHY    : for split physics (one part at t-Dt, one part at t+Dt);
!             can be used only at ECMWF with the ECMWF package.

! ----- way of computing initial state:
! N2DINI    : 1 initialization for 2D with Haurwitz wave
!           : 2 initialization for 2D with real fields
! N3DINI    : 1 initialization for 3D with standard atmosphere
!             (not available since cycle 12)
!           : 0 initialization for 3D with real fields

! ----- control of variables which are transformed in spectral space:
! LSPRT     : .T.: if R*T/Rd "virtual temperature" as spectral variable

! ----- diagnostics and frequencies:
! NSPPR     : 0: no spectrum printed in spnorm; only global norms averaged
!             on the vertical are printed
!           : 1: no spectrum printed in spnorm; only global norms averaged
!             on the vertical and global norms on each layer are printed
!           : 2: both total wavenumber spectrum and zonal wavenumber spectrum
!             are printed
! NFRPOS    : frequency of post-processing events (time-steps)
! NFRISP    : frequency of isp (animation !) events (time-steps)
! NFRCO     : frequency of coupled fields (time-steps)
! NFRCORM   : mass correction frequency (>0 time-steps, <0 hours, 0 no
!             correction)
! NFRHIS    : frequency of history write_ups (time-steps)
! NFRSFXHIS   : frequency of history write_ups (time-steps) for surface
! NFRMASSCON: frequency of mass conservation fixup (time-steps)
! N6BINS    : maximum frequency of 6 hour min/max post processing events (hours)  
! NFRGDI    : frequency of grid-point space diagnostics
! NFRSDI    : frequency of spectral space diagnostics
! NFRDHFG   : write-up frequency of global DDH
! NFRDHFZ   : write-up frequency of zonal DDH
! NFRDHFD   : write-up frequency of "limited domain" DDH
! NFRDHP    : write-up frequency of DDH files
! NPOSTS    : array containing postprocessing steps
! NPOSTSMIN : array containing postprocessing steps (minutes)
! NPISPS    : array containing isp (animation !) steps
! NHISTS    : array containing history write-up steps
! NHISTSMIN : array containing history write-up steps (minutes)
! NSFXHISTS   : array containing history write-up steps for surface
! NSFXHISTSMIN : array containing history write-up steps (minutes) for surface
! NMASSCONS : array containing mass conservation fixup steps
! NGDITS    : array containing grid point diagnostics steps
! NSDITS    : array containing spectral diagnostics steps
! NDHFGTS   : array containing write out steps for global DDH
! NDHFZTS   : array containing write out steps for zonal means DDH
! NDHFDTS   : array containing write out steps for limited areas DDH
! NDHPTS    : array containing write out steps for DDH
! Explanation for N[XXX]TS:
!             1) if N[XXX]TS(0)=0 action if MOD(JSTEP,NFR[XXX])=0
!             2) if N[XXX]TS(0)>0 N[XXX]TS(0) significant numbers in 
!                N[XXX]TS are then considered and:
!                action for JSTEP=N[XXX]TS(.)*NFR[XXX]
!             3) IF N[XXX]TS(0)<0
!                action for JSTEP=(N[XXX]TS(.)/DELTAT)*NFR[XXX]

!========== ASSIMILATION SWITCHES =============================================

! L4DVAR  : .T. => 4DVAR (occurs only if NSTOP>1 and NCONF=131)
! LNOBGON : .F. if eq. .T. no stop if unsuccessful reading of CMA files
! LCANARI : .T. = term to control French OI
! NTASKS_CANARI : Number of tasks for CANARI.
! LCASIG  : .T. = interpolation of the model errors (French OI)
! LGUESS  : .T. = term of first guess included
! LOBS    : .T. = term of observations included
! LSIMOB  : .T. = if simulated observations
! LOBSC1  : .T. = term of observations included in configuration 1
! LSCREEN : .T. = observation screening for variational assimilation
! LSCREEN_OPENMP : .T. = 4DVAR screening runs in OpenMP-parallel mode over timeslots
! L_SPLIT_SCREEN .T. = to split screenng
! L_SCREEN_CALL: .T. = call to screening routine SCREEN
! LMONITORING : .T. = monitor data which are not in 4dvar
! LOBSREF : .T. = comparison to observation for the trajectory (NCONF=131)
! LIFSMIN : .T. = if running minimisation
! LIFSTRAJ: .T. = if running high resolution trajectory integration
! NCNTVAR : Definition of the control variable of a variational job.
!           = 1 ===> control variables are model variables
!           = 2 ===> control variables are normalized departures of
!                    model variables from the background field
!           = 3 ===> ..........
! LOLDPP  : .T. use "old" p.p. of T,Q,U and V
! NSTEPINI: Initial step in hours for the initial conditions
!           at the beginning of 4D-Var trajectory (usually 3 hours).
!           It is used to update the step while saving the FCs along
!           the first trajectory.
! NINTERPTRAJ : Interpolation method applied to increments
! NINTERPINCR : Interpolation method applied to increments
!           = 1 ===> Bi-linear interpolation (default)
!           = 2 ===> Bi-cubic Rinterpolation
!           = 3 ===> Conserving interpolation
! For Conserving interpolation the style is defined by:
!
!     LINFLAT       .T. Inflation of perturbed backgrounds (6h fc wrt beg. of window)
!     LINFLP9       .T. Inflation of 9h fcsts (wrt beg. of window)
!     LINFLF1       .T. Neutral inflation (i.e. inflation factor =1)

!========== ECMWF / METEO-FRANCE SWITCHES ======================================

!*  Control variables for the default set-up DO NOT USE THEM OUT OF THE SET-UP !ub

! LECMWF  : .T. = set-up default values used at ECMWF
!           .F. = set-up default values used at METEO-FRANCE

!========== ALADIN SWITCHES ===================================================

! LELAM   : .T. = limited area model with coupling or fully biperiodic
!                 model (torus)
!           .F. = global model (sphere)
! LRPLANE : .T. = plane cartesian geometry
!           .F. = spherical geometry

!========== AROME SWITCH ======================================================

! LAROME  : .T. = AROME limited area model

!========== SURFEX SWITCH ======================================================

! LCALLSFX : .F. = To prevent SURFEX to be called twice at NSTEP=0
! LSFXLSM  : .F. = To prevent to use the LandSea mask from SURFEX      

!========== ECMWF Single Column Model =========================================
! LSFCFLX : .T. = forcing with surface fluxes (latent and sensible).
! REXTSHF : externally supplied sensible heat flux [W/m^2]
! REXTLHF : externally supplied latent   heat flux [W/m^2]
! LROUGH  : .T. = surface roughness length is externally specified 
! REXTZ0M : externally supplied roughness length for momentum [m]
! REXTZ0H : externally supplied roughness length for heat [m]

!========== FORCING SWITCH ====================================================
! LSFORC - switch to activate the large scale forcings in setup and cpg
! LSFORCS - switch to activate the surface forcings
!==============================================================================

!========== Coupled Runs with OASIS4 ==========================================
! LCOUPLO4 -  coupled runs, i.e. get communicator from prism  
! LCOUPLO4_ENV -  coupled runs, environment variable  
!----------------------------------------------

! * Parameters:
INTEGER(KIND=JPIM), PARAMETER :: JPNPST=960

! * Technical switches:
INTEGER(KIND=JPIM) :: NCONF
INTEGER(KIND=JPIM) :: NFPOS
LOGICAL :: LR3D
LOGICAL :: LR2D
LOGICAL :: LRSHW
LOGICAL :: LRVEQ
LOGICAL :: LALLOPR
LOGICAL :: LFDBOP
LOGICAL :: LECFPOS
LOGICAL :: LGRBOP
LOGICAL :: LFBDAP
LOGICAL :: LARPEGEF
LOGICAL :: LARPEGEF_TRAJHR
LOGICAL :: LARPEGEF_TRAJBG
LOGICAL :: LARPEGEF_RDGP_INIT
LOGICAL :: LARPEGEF_RDGP_TRAJHR
LOGICAL :: LARPEGEF_RDGP_TRAJBG
CHARACTER (LEN = 120) ::  CNDISPP
CHARACTER (LEN = 128) ::  CFPNCF
LOGICAL :: LFPART2
CHARACTER (LEN = 120) ::  CFDIRLST
CHARACTER (LEN = 120) ::  CNPPATH
LOGICAL :: LRETCFOU
LOGICAL :: LWRTCFOU
LOGICAL :: LRFOUTCNORM(4)
LOGICAL :: LRGPTCNORM(3)
LOGICAL :: LNF
LOGICAL :: LSMSSIG
CHARACTER (LEN = 25) :: CEVENT
CHARACTER (LEN = 25) :: CMETER
LOGICAL :: LOPDIS
INTEGER(KIND=JPIM) :: NCYCLE=145
CHARACTER (LEN = 16) ::  CNMEXP
CHARACTER (LEN = 2) ::  CFCLASS
CHARACTER (LEN = 2) ::  CTYPE
LOGICAL :: LBACKG
LOGICAL :: LSPBSBAL
LOGICAL :: LMINIM
LOGICAL :: LIOLEVG
INTEGER(KIND=JPIM) :: NUNDEFLD
LOGICAL :: LREFOUT
LOGICAL :: LREFGEN

! * Model switches:
LOGICAL :: LSLAG
LOGICAL :: LTWOTL
LOGICAL :: LREGETA
LOGICAL :: LVFE_REGETA
INTEGER(KIND=JPIM) :: NQUAD
LOGICAL :: LNHDYN
LOGICAL :: LSLPHY
INTEGER(KIND=JPIM) :: N2DINI
INTEGER(KIND=JPIM) :: N3DINI
LOGICAL :: LSPRT
INTEGER(KIND=JPIM) :: NSPPR
INTEGER(KIND=JPIM) :: NFRPOS
INTEGER(KIND=JPIM) :: NFRISP
INTEGER(KIND=JPIM) :: NFRCORM
INTEGER(KIND=JPIM) :: NFRCO
INTEGER(KIND=JPIM) :: NFRHIS
INTEGER(KIND=JPIM) :: NFRSFXHIS
INTEGER(KIND=JPIM) :: NFRMASSCON
INTEGER(KIND=JPIM) :: NFRGDI
INTEGER(KIND=JPIM) :: NFRSDI
INTEGER(KIND=JPIM) :: NFRDHFG
INTEGER(KIND=JPIM) :: NFRDHFZ
INTEGER(KIND=JPIM) :: NFRDHFD
INTEGER(KIND=JPIM) :: NFRDHP
INTEGER(KIND=JPIM) :: N6BINS
INTEGER(KIND=JPIM) :: NPOSTS(0:JPNPST)
INTEGER(KIND=JPIM) :: NPOSTSMIN(0:JPNPST)
INTEGER(KIND=JPIM) :: NPISPS(0:JPNPST)
INTEGER(KIND=JPIM) :: NHISTS(0:JPNPST)
INTEGER(KIND=JPIM) :: NHISTSMIN(0:JPNPST)
INTEGER(KIND=JPIM) :: NSFXHISTS(0:JPNPST)
INTEGER(KIND=JPIM) :: NSFXHISTSMIN(0:JPNPST)
INTEGER(KIND=JPIM) :: NMASSCONS(0:JPNPST)
INTEGER(KIND=JPIM) :: NGDITS(0:JPNPST)
INTEGER(KIND=JPIM) :: NSDITS(0:JPNPST)
INTEGER(KIND=JPIM) :: NDHFGTS(0:JPNPST)
INTEGER(KIND=JPIM) :: NDHFZTS(0:JPNPST)
INTEGER(KIND=JPIM) :: NDHFDTS(0:JPNPST)
INTEGER(KIND=JPIM) :: NDHPTS(0:JPNPST)

! * Assimilation:
LOGICAL :: L4DVAR=.FALSE.
LOGICAL :: LNOBGON
LOGICAL :: LCANARI
INTEGER(KIND=JPIM) :: NTASKS_CANARI
LOGICAL :: LCASIG
LOGICAL :: LGUESS
LOGICAL :: LOBS
LOGICAL :: LSIMOB
LOGICAL :: LOBSC1
LOGICAL :: LSCREEN
LOGICAL :: LSCREEN_OPENMP
LOGICAL :: L_SPLIT_SCREEN
LOGICAL :: L_SCREEN_CALL
LOGICAL :: LMONITORING
LOGICAL :: LOBSREF
LOGICAL :: LIFSMIN
LOGICAL :: LIFSTRAJ
INTEGER(KIND=JPIM) :: NCNTVAR
LOGICAL :: LOLDPP
INTEGER(KIND=JPIM) :: NSTEPINI
INTEGER(KIND=JPIM) :: NINTERPTRAJ
INTEGER(KIND=JPIM) :: NINTERPINCR
LOGICAL :: LINFLAT
LOGICAL :: LINFLP9
LOGICAL :: LINFLF1

! * ECMWF / METEO-FRANCE:
LOGICAL :: LECMWF

! * ALADIN:
LOGICAL :: LELAM
LOGICAL :: LRPLANE

! * AROME:
LOGICAL :: LAROME

! * SURFEX:
LOGICAL :: LCALLSFX
LOGICAL :: LSFXLSM  

! * ECMWF Single Column Model:
LOGICAL :: LSFCFLX
REAL(KIND=JPRB) :: REXTSHF
REAL(KIND=JPRB) :: REXTLHF
LOGICAL :: LROUGH
REAL(KIND=JPRB) :: REXTZ0M
REAL(KIND=JPRB) :: REXTZ0H

! * coupled runs OASIS4
LOGICAL :: LCOUPLO4
LOGICAL :: LCOUPLO4_ENV

! * FORCING
LOGICAL :: LSFORC
LOGICAL :: LSFORCS

! * SPECTRAL/GRID-POINT IO

! Write spectral data in GP using WRSPECA_GP
LOGICAL :: LWRSPECA_GP = .FALSE.
! Read spectral data in GP using SUSPECA_GP
LOGICAL :: LSUSPECA_GP = .FALSE.
! Write spectral data in U/V GP representation
LOGICAL :: LWRSPECA_GP_UV = .FALSE.
! Read spectral data in U/V GP representation
LOGICAL :: LSUSPECA_GP_UV = .FALSE.

!     ------------------------------------------------------------------
END MODULE YOMCT0
