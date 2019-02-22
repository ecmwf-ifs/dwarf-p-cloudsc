#include "cloudsc_driver.h"

void cloudsc_driver(int numthreads, int numcols, int nproma) {

    double *tend_tmp_u;
    double *tend_tmp_v;
    double *tend_tmp_t;
    double *tend_tmp_q;
    double *tend_tmp_o3;
    double *tend_tmp_a;
    double *tend_tmp_cld;

    double *tend_loc_u;
    double *tend_loc_v;
    double *tend_loc_t;
    double *tend_loc_q;
    double *tend_loc_o3;
    double *tend_loc_a;
    double *tend_loc_cld;

    double *tend_cml_u;
    double *tend_cml_v;
    double *tend_cml_t;
    double *tend_cml_q;
    double *tend_cml_o3;
    double *tend_cml_a;
    double *tend_cml_cld;

    double ptsphy;            //! Physics timestep

    double *plcrit_aer, *zplcrit_aer;
    double *picrit_aer, *zpicrit_aer;
    double *pre_ice, *zpre_ice;
    double *pccn, *zpccn;       //! liquid cloud condensation nuclei
    double *pnice, *zpnice;     //! ice number concentration (cf. CCN)
    double *pt, *zpt;           //! T at start of callpar
    double *pq, *zpq;           //! Q at start of callpar
    double *pvfa, *zpvfa;       //! CC from VDF scheme
    double *pvfl, *zpvfl;       //! Liq from VDF scheme
    double *pvfi, *zpvfi;       //! Ice from VDF scheme
    double *pdyna, *zpdyna;     //! CC from Dynamics
    double *pdynl, *zpdynl;     //! Liq from Dynamics
    double *pdyni, *zpdyni;     //! Liq from Dynamics
    double *phrsw, *zphrsw;     //! Short-wave heating rate
    double *phrlw, *zphrlw;     //! Long-wave heating rate
    double *pvervel, *zpvervel; //! Vertical velocity
    double *pap, *zpap;         //! Pressure on full levels
    double *paph, *zpaph;       //! Pressure on half levels
    double *plsm, *zplsm;           //! Land fraction (0-1) 
    int    *ldcum, *llcum;          //! Convection active
    int    *ktype, *itype;          //! Convection type 0,1,2
    double *plu, *zplu;         //! Conv. condensate
    double *plude, *zplude;     //! Conv. detrained water 
    double *plude_tmp;
    double *psnde, *zpsnde;     //! Conv. detrained snow
    double *pmfu, *zpmfu;       //! Conv. mass flux up
    double *pmfd, *zpmfd;       //! Conv. mass flux down
    double *pa, *zpa;           //! Original Cloud fraction (t)

    double *pclv, zpclv;
    double *psupsat, zpsupsat;

    double *pcovptot, *zpcovptot;   //! Precip fraction
    double *prainfrac_toprfz, *zprainfrac_toprfz;
    double *pfsqlf, *zpfsqlf;       //! Flux of liquid
    double *pfsqif, *zpfsqif;       //! Flux of ice
    double *pfcqlng, *zpfcqlng;     //! -ve corr for liq
    double *pfcqnng, *zpfcqnng;     //! -ve corr for ice
    double *pfsqrf, *zpfsqrf;       //! Flux diagnostics
    double *pfsqsf, *zpfsqsf;       //!    for DDH, generic
    double *pfcqrng, *zpfcqrng;     //! rain
    double *pfcqsng, *zpfcqsng;     //! snow
    double *pfsqltur, *zpfsqltur;   //! liquid flux due to VDF
    double *pfsqitur, *zpfsqitur;   //! ice flux due to VDF
    double *pfplsl, *zpfplsl;       //! liq+rain sedim flux
    double *pfplsn, *zpfplsn;       //! ice+snow sedim flux
    double *pfhpsl, *zpfhpsl;       //! Enthalpy flux for liq
    double *pfhpsn, *zpfhpsn;       //! Enthalpy flux for ice

    /* Define or query data dimensions from input file */
    int klon, nlev;
    int kidia, kfdia;
    int jkglo,ibl,icend;

    nclv = 5;      // number of microphysics variables
    ncldql = 1;    // liquid cloud water
    ncldqi = 2;    // ice cloud water
    ncldqr = 3;    // rain water
    ncldqs = 4;    // snow
    ncldqv = 5;    // vapour

    yrecldp = malloc(sizeof(struct TECLDP));

    query_state(&klon, &nlev);

    tend_loc_u   = (double*) malloc( sizeof(double) * numcols*nlev ); 
    tend_loc_v   = (double*) malloc( sizeof(double) * numcols*nlev );
    tend_loc_t   = (double*) malloc( sizeof(double) * numcols*nlev );
    tend_loc_o3  = (double*) malloc( sizeof(double) * numcols*nlev );
    tend_loc_q   = (double*) malloc( sizeof(double) * numcols*nlev );
    tend_loc_a   = (double*) malloc( sizeof(double) * numcols*nlev );
    tend_loc_cld = (double*) malloc( sizeof(double) * numcols*nlev*nclv ); 

    tend_cml_u   = (double*) malloc( sizeof(double) * numcols*nlev ); 
    tend_cml_v   = (double*) malloc( sizeof(double) * numcols*nlev ); 
    tend_cml_t   = (double*) malloc( sizeof(double) * numcols*nlev ); 
    tend_cml_o3  = (double*) malloc( sizeof(double) * numcols*nlev );
    tend_cml_q   = (double*) malloc( sizeof(double) * numcols*nlev );
    tend_cml_a   = (double*) malloc( sizeof(double) * numcols*nlev );
    tend_cml_cld = (double*) malloc( sizeof(double) * numcols*nlev*nclv );

    tend_tmp_u   = (double*) malloc( sizeof(double) * numcols*nlev );
    tend_tmp_v   = (double*) malloc( sizeof(double) * numcols*nlev );
    tend_tmp_t   = (double*) malloc( sizeof(double) * numcols*nlev );
    tend_tmp_o3  = (double*) malloc( sizeof(double) * numcols*nlev );
    tend_tmp_q   = (double*) malloc( sizeof(double) * numcols*nlev );
    tend_tmp_a   = (double*) malloc( sizeof(double) * numcols*nlev );
    tend_tmp_cld = (double*) malloc( sizeof(double) * numcols*nlev*nclv );


    plcrit_aer = (double*) malloc( sizeof(double) * numcols*nlev );
    picrit_aer = (double*) malloc( sizeof(double) * numcols*nlev );
    pre_ice    = (double*) malloc( sizeof(double) * numcols*nlev );
    pccn       = (double*) malloc( sizeof(double) * numcols*nlev );
    pnice      = (double*) malloc( sizeof(double) * numcols*nlev );
    pt         = (double*) malloc( sizeof(double) * numcols*nlev );
    pq         = (double*) malloc( sizeof(double) * numcols*nlev );
    pvfa       = (double*) malloc( sizeof(double) * numcols*nlev );
    pvfl       = (double*) malloc( sizeof(double) * numcols*nlev );
    pvfi       = (double*) malloc( sizeof(double) * numcols*nlev );
    pdyna      = (double*) malloc( sizeof(double) * numcols*nlev );
    pdynl      = (double*) malloc( sizeof(double) * numcols*nlev );
    pdyni      = (double*) malloc( sizeof(double) * numcols*nlev );
    phrsw      = (double*) malloc( sizeof(double) * numcols*nlev );
    phrlw      = (double*) malloc( sizeof(double) * numcols*nlev );
    pvervel    = (double*) malloc( sizeof(double) * numcols*nlev );
    pap        = (double*) malloc( sizeof(double) * numcols*nlev );
    paph       = (double*) malloc( sizeof(double) * numcols*(nlev+1) );
    plsm       = (double*) malloc( sizeof(double) * numcols );
    ldcum      = (int*) malloc( sizeof(int) * numcols );
    ktype      = (int*) malloc( sizeof(int) * numcols );
    plu        = (double*) malloc( sizeof(double) * numcols*nlev );
    plude      = (double*) malloc( sizeof(double) * numcols*nlev );
    psnde      = (double*) malloc( sizeof(double) * numcols*nlev );
    pmfu       = (double*) malloc( sizeof(double) * numcols*nlev );
    pmfd       = (double*) malloc( sizeof(double) * numcols*nlev );
    pa         = (double*) malloc( sizeof(double) * numcols*nlev );
    pclv       = (double*) malloc( sizeof(double) * numcols*nlev*nclv );
    psupsat    = (double*) malloc( sizeof(double) * numcols*nlev );
    pcovptot   = (double*) malloc( sizeof(double) * numcols*nlev );


    prainfrac_toprfz = (double*) malloc( sizeof(double) * numcols );
    pfsqlf    = (double*) malloc( sizeof(double) * numcols*(nlev+1) );
    pfsqif    = (double*) malloc( sizeof(double) * numcols*(nlev+1) );
    pfcqnng   = (double*) malloc( sizeof(double) * numcols*(nlev+1) );
    pfcqlng   = (double*) malloc( sizeof(double) * numcols*(nlev+1) );
    pfsqrf    = (double*) malloc( sizeof(double) * numcols*(nlev+1) );
    pfsqsf    = (double*) malloc( sizeof(double) * numcols*(nlev+1) );
    pfcqrng   = (double*) malloc( sizeof(double) * numcols*(nlev+1) );
    pfcqsng   = (double*) malloc( sizeof(double) * numcols*(nlev+1) );
    pfsqltur  = (double*) malloc( sizeof(double) * numcols*(nlev+1) );
    pfsqitur  = (double*) malloc( sizeof(double) * numcols*(nlev+1) );
    pfplsl    = (double*) malloc( sizeof(double) * numcols*(nlev+1) );
    pfplsn    = (double*) malloc( sizeof(double) * numcols*(nlev+1) );
    pfhpsl    = (double*) malloc( sizeof(double) * numcols*(nlev+1) );
    pfhpsn    = (double*) malloc( sizeof(double) * numcols*(nlev+1) );

    load_state(klon, nlev, nclv, numcols, nproma, &ptsphy, plcrit_aer, picrit_aer,
	       pre_ice, pccn, pnice, pt, pq,
	       tend_cml_t, tend_cml_q, tend_cml_a, tend_cml_cld,
	       tend_tmp_t, tend_tmp_q, tend_tmp_a, tend_tmp_cld,
	       pvfa, pvfl, pvfi, pdyna, pdynl, pdyni,
	       phrsw, phrlw, pvervel, pap, paph, plsm, ldcum, ktype, plu,
	       plude, psnde, pmfu, pmfd, pa, pclv, psupsat);

    cloudsc_c(1, numcols, numcols, nlev, ptsphy,  pt, pq,  
              tend_cml_t, tend_cml_q, tend_cml_a, tend_cml_cld,
              tend_tmp_t, tend_tmp_q, tend_tmp_a, tend_tmp_cld,
              tend_loc_t, tend_loc_q, tend_loc_a, tend_loc_cld,
              pvfa,  pvfl,  pvfi, 
              pdyna,  pdynl,  pdyni,  
              phrsw, phrlw,  pvervel, 
              pap,  paph,  plsm, ldcum, ktype, 
              plu, plude, psnde,  pmfu,  pmfd,  
              pa, pclv, psupsat,  
              plcrit_aer,  picrit_aer,  pre_ice, pccn,  pnice,
              pcovptot,  prainfrac_toprfz,  pfsqlf,  
              pfsqif, pfcqnng,  pfcqlng,
              pfsqrf,  pfsqsf,  pfcqrng,  
              pfcqsng, pfsqltur, pfsqitur,
              pfplsl,  pfplsn,  pfhpsl,  pfhpsn);

    printf("finished cloudsc, now checking results\n");

    /* check_results(numcols, nlev, nclv,    plude ,  pcovptot ,  prainfrac_toprfz ,  pfsqlf ,  pfsqif , */
    /*           pfcqlng ,  pfcqnng ,  pfsqrf ,  pfsqsf ,  pfcqrng ,  pfcqsng , */
    /*           pfsqltur ,  pfsqitur ,  pfplsl ,  pfplsn ,  pfhpsl ,  pfhpsn , */
    /*           tend_loc_a ,  tend_loc_q ,  tend_loc_t ,  tend_loc_cld ) ; */

    free(plcrit_aer); // ALLOCATE(PLCRIT_AER(KLON,KLEV)) 
    free(picrit_aer); // ALLOCATE(PICRIT_AER(KLON,KLEV)) 
    free(pre_ice);    // ALLOCATE(PRE_ICE(KLON,KLEV)) 
    free(pccn);
    free(pnice);
    free(pt);
    free(pq);
    free(pvfa);
    free(pvfl);
    free(pvfi);
    free(pdyna);
    free(pdynl);
    free(pdyni);
    free(phrsw);
    free(phrlw);
    free(pvervel);
    free(pap);
    free(paph); // ALLOCATE(PAPH(KLON,KLEV+1))
    free(plsm);
    free(ldcum);
    free(ktype);
    free(plu);
    free(plude);
    free(psnde);
    free(pmfu);
    free(pmfd);
    free(pa);
    free(pclv);
    free(psupsat);
    free(pcovptot);
    free(tend_loc_u);
    free(tend_loc_v);
    free(tend_loc_t);
    free(tend_loc_o3);
    free(tend_loc_q);
    free(tend_loc_a);
    free(tend_loc_cld);
    free(tend_tmp_u);
    free(tend_tmp_v);
    free(tend_tmp_t);
    free(tend_tmp_o3);
    free(tend_tmp_q);
    free(tend_tmp_a);
    free(tend_tmp_cld);
    free(tend_cml_u);
    free(tend_cml_v);
    free(tend_cml_t);
    free(tend_cml_o3);
    free(tend_cml_q);
    free(tend_cml_a);
    free(tend_cml_cld);

    free(prainfrac_toprfz);
    free(pfsqlf);
    free(pfsqif);
    free(pfcqnng);
    free(pfcqlng);
    free(pfsqrf);
    free(pfsqsf);
    free(pfcqrng);
    free(pfcqsng);
    free(pfsqltur);
    free(pfsqitur);
    free(pfplsl);
    free(pfplsn);
    free(pfhpsl);
    free(pfhpsn);

    free(yrecldp);
}
