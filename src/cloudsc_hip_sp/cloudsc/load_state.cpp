/*
 * (C) Copyright 1988- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "load_state.h"
//#include "yomcst_c.hpp"
#include <iostream>

#include <math.h>
#include "serialbox-c/Serialbox.h"

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

/* Query sizes and dimensions of state arrays */
void query_state(int *klon, int *klev)
{
  serialboxSerializer_t* serializer = serialboxSerializerCreate(Read, "./data", "input", "Binary");
  serialboxMetainfo_t* globalMetainfo = serialboxSerializerGetGlobalMetainfo(serializer);

  *klon = serialboxMetainfoGetInt32(globalMetainfo, "KLON");
  *klev = serialboxMetainfoGetInt32(globalMetainfo, "KLEV");

  serialboxMetainfoDestroy(globalMetainfo);
  serialboxSerializerDestroy(serializer);
}

void expand_1d(float *buffer, float *field_in, int nlon, int nproma, int ngptot, int nblocks)
{
  int b, i, buf_start_idx, buf_idx;

#pragma omp parallel for default(shared) private(b, i, buf_start_idx, buf_idx)
  for (b = 0; b < nblocks; b++) {
    buf_start_idx = ((b)*nproma) % nlon;
    for (i = 0; i < nproma; i++) {
      buf_idx = (buf_start_idx + i) % nlon;
      field_in[b*nproma+i] = buffer[buf_idx];
    }
  }
}


void expand_1d_int(int *buffer, int *field_in, int nlon, int nproma, int ngptot, int nblocks)
{
  int b, i, buf_start_idx, buf_idx;

  #pragma omp parallel for default(shared) private(b, i, buf_start_idx, buf_idx)
  for (b = 0; b < nblocks; b++) {
    buf_start_idx = ((b)*nproma) % nlon;
    for (i = 0; i < nproma; i++) {
      buf_idx = (buf_start_idx + i) % nlon;
      field_in[b*nproma+i] = buffer[buf_idx];
    }
  }
}


void expand_2d(float *buffer_in, float *field_in, int nlon, int nlev, int nproma, int ngptot, int nblocks)
{
  int b, l, i, buf_start_idx, buf_idx;

  #pragma omp parallel for default(shared) private(b, buf_start_idx, buf_idx, l, i)
  for (b = 0; b < nblocks; b++) {
    buf_start_idx = ((b)*nproma) % nlon;
    for (i = 0; i < nproma; i++) {
      for (l = 0; l < nlev; l++) {
        buf_idx = (buf_start_idx + i) % nlon;
        field_in[b*nlev*nproma+l*nproma+i] = buffer_in[l*nlon+buf_idx];
      }
    }
  }
}

void expand_3d(float *buffer_in, float *field_in, int nlon, int nlev, int nclv, int nproma, int ngptot, int nblocks)
{
  int b, l, c, i, buf_start_idx, buf_idx;

#pragma omp parallel for default(shared) private(b, buf_start_idx, buf_idx, l, i)
  for (b = 0; b < nblocks; b++) {
    buf_start_idx = ((b)*nproma) % nlon;
    for (i = 0; i < nproma; i++) {
      for (c = 0; c < nclv; c++) {
        for (l = 0; l < nlev; l++) {
          buf_idx = (buf_start_idx + i) % nlon;
          field_in[b*nclv*nlev*nproma+c*nlev*nproma+l*nproma+i] = buffer_in[c*nlev*nlon+l*nlon+buf_idx];
        }
      }
    }
  }
}


void load_and_expand_1d(serialboxSerializer_t *serializer, serialboxSavepoint_t* savepoint,
    const char *name, int nlon, int nproma, int ngptot, int nblocks, float *field)
{
  double dbl_buffer[nlon];
  float buffer[nlon];
  int strides[1] = {1};

  serialboxSerializerRead(serializer, name, savepoint, dbl_buffer, strides, 1);
  for (int i=0;i<nlon;i++) {
    buffer[i] = (float) dbl_buffer[i];
  }
  expand_1d((float *)buffer, field, nlon, nproma, ngptot, nblocks);
}

void load_and_expand_1d_int(serialboxSerializer_t *serializer, serialboxSavepoint_t* savepoint,
    const char *name, int nlon, int nproma, int ngptot, int nblocks, int *field)
{
  int buffer[nlon]; // TODO: double?
  int strides[1] = {1};

  serialboxSerializerRead(serializer, name, savepoint, buffer, strides, 1);
  expand_1d_int((int *)buffer, field, nlon, nproma, ngptot, nblocks);
}

void load_and_expand_2d(serialboxSerializer_t *serializer, serialboxSavepoint_t* savepoint,
    const char *name, int nlon, int nlev, int nproma, int ngptot, int nblocks, float *field)
{
  double dbl_buffer[nlev][nlon];
  float buffer[nlev][nlon];
  int strides[2] = {1, nlon};

  serialboxSerializerRead(serializer, name, savepoint, dbl_buffer, strides, 2);
  for (int j=0;j<nlev;j++) {
  for (int i=0;i<nlon;i++) {
    buffer[j][i] = (float) dbl_buffer[j][i];
  }
  }
  expand_2d((float *)buffer, field, nlon, nlev, nproma, ngptot, nblocks);
}

void load_and_expand_3d(serialboxSerializer_t *serializer, serialboxSavepoint_t* savepoint,
    const char *name, int nlon, int nlev, int nclv, int nproma, int ngptot, int nblocks, float *field)
{
  double dbl_buffer[nclv][nlev][nlon];
  float buffer[nclv][nlev][nlon];
  int strides[3] = {1, nlon, nlev*nlon};

  serialboxSerializerRead(serializer, name, savepoint, dbl_buffer, strides, 3);
  for (int k=0;k<nclv;k++) {
  for (int j=0;j<nlev;j++) {
  for (int i=0;i<nlon;i++) {
    buffer[k][j][i] = (float) dbl_buffer[k][j][i];
  }
  }
  }
  expand_3d((float *)buffer, field, nlon, nlev, nclv, nproma, ngptot, nblocks);
}


/* Read input state into memory */
void load_state(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma, 
    float* ptsphy, float* plcrit_aer, float* picrit_aer,
    float* pre_ice, float* pccn, float* pnice, float* pt, float* pq, 
    float* tend_cml_t, float* tend_cml_q, float* tend_cml_a, float* tend_cml_cld,
    float* tend_tmp_t, float* tend_tmp_q, float* tend_tmp_a, float* tend_tmp_cld,
    float* pvfa, float* pvfl, float* pvfi, float* pdyna, float* pdynl, float* pdyni, 
    float* phrsw, float* phrlw, float* pvervel, float* pap, float* paph, float* plsm,
    int* ktype, float* plu, float* plude, float* psnde, float* pmfu,
    float* pmfd, float* pa, float* pclv, float* psupsat, struct TECLDP* yrecldp,
    float* rg, float* rd, float* rcpd, float* retv, float* rlvtt, float* rlstt, 
    float* rlmlt, float* rtt, float* rv, float* r2es, float* r3les, float* r3ies,
    float* r4les, float* r4ies, float* r5les, float* r5ies, float* r5alvcp, float* r5alscp,
    float* ralvdcp, float* ralsdcp, float* ralfdcp, float* rtwat, 
    float* rtice, float* rticecu, float* rtwat_rtice_r, float *rtwat_rticecu_r,
    float* rkoop1, float* rkoop2 )
{
  serialboxSerializer_t* serializer = serialboxSerializerCreate(Read, "./data", "input", "Binary");
  serialboxMetainfo_t* metainfo = serialboxSerializerGetGlobalMetainfo(serializer);
  serialboxSavepoint_t** savepoints = serialboxSerializerGetSavepointVector(serializer);
  serialboxSavepoint_t* savepoint = savepoints[0];

  int nblocks = (ngptot / nproma) + min(ngptot % nproma, 1);

  load_and_expand_2d(serializer, savepoint, "PLCRIT_AER", nlon, nlev, nproma, ngptot, nblocks, plcrit_aer);
  load_and_expand_2d(serializer, savepoint, "PICRIT_AER", nlon, nlev, nproma, ngptot, nblocks, picrit_aer);
  load_and_expand_2d(serializer, savepoint, "PRE_ICE", nlon, nlev, nproma, ngptot, nblocks, pre_ice);
  load_and_expand_2d(serializer, savepoint, "PCCN", nlon, nlev, nproma, ngptot, nblocks, pccn);
  load_and_expand_2d(serializer, savepoint, "PNICE", nlon, nlev, nproma, ngptot, nblocks, pnice);
  load_and_expand_2d(serializer, savepoint, "PT", nlon, nlev, nproma, ngptot, nblocks, pt);
  load_and_expand_2d(serializer, savepoint, "PQ", nlon, nlev, nproma, ngptot, nblocks, pq);
  load_and_expand_2d(serializer, savepoint, "TENDENCY_CML_T", nlon, nlev, nproma, ngptot, nblocks, tend_cml_t);
  load_and_expand_2d(serializer, savepoint, "TENDENCY_CML_Q", nlon, nlev, nproma, ngptot, nblocks, tend_cml_q);
  load_and_expand_2d(serializer, savepoint, "TENDENCY_CML_A", nlon, nlev, nproma, ngptot, nblocks, tend_cml_a);
  load_and_expand_3d(serializer, savepoint, "TENDENCY_CML_CLD", nlon, nlev, nclv, nproma, ngptot, nblocks, tend_cml_cld);
  load_and_expand_2d(serializer, savepoint, "TENDENCY_TMP_T", nlon, nlev, nproma, ngptot, nblocks, tend_tmp_t);
  load_and_expand_2d(serializer, savepoint, "TENDENCY_TMP_Q", nlon, nlev, nproma, ngptot, nblocks, tend_tmp_q);
  load_and_expand_2d(serializer, savepoint, "TENDENCY_TMP_A", nlon, nlev, nproma, ngptot, nblocks, tend_tmp_a);
  load_and_expand_3d(serializer, savepoint, "TENDENCY_TMP_CLD", nlon, nlev, nclv, nproma, ngptot, nblocks, tend_tmp_cld);
  load_and_expand_2d(serializer, savepoint, "PVFA", nlon, nlev, nproma, ngptot, nblocks, pvfa);
  load_and_expand_2d(serializer, savepoint, "PVFL", nlon, nlev, nproma, ngptot, nblocks, pvfl);
  load_and_expand_2d(serializer, savepoint, "PVFI", nlon, nlev, nproma, ngptot, nblocks, pvfi);
  load_and_expand_2d(serializer, savepoint, "PDYNA", nlon, nlev, nproma, ngptot, nblocks, pdyna);
  load_and_expand_2d(serializer, savepoint, "PDYNL", nlon, nlev, nproma, ngptot, nblocks, pdynl);
  load_and_expand_2d(serializer, savepoint, "PDYNI", nlon, nlev, nproma, ngptot, nblocks, pdyni);
  load_and_expand_2d(serializer, savepoint, "PHRSW", nlon, nlev, nproma, ngptot, nblocks, phrsw);
  load_and_expand_2d(serializer, savepoint, "PHRLW", nlon, nlev, nproma, ngptot, nblocks, phrlw);
  load_and_expand_2d(serializer, savepoint, "PVERVEL", nlon, nlev, nproma, ngptot, nblocks, pvervel);
  load_and_expand_2d(serializer, savepoint, "PAP", nlon, nlev, nproma, ngptot, nblocks, pap);
  load_and_expand_2d(serializer, savepoint, "PAPH", nlon, nlev+1, nproma, ngptot, nblocks, paph);
  load_and_expand_1d(serializer, savepoint, "PLSM", nlon, nproma, ngptot, nblocks, plsm);
  load_and_expand_1d_int(serializer, savepoint, "KTYPE", nlon, nproma, ngptot, nblocks, ktype);
  load_and_expand_2d(serializer, savepoint, "PLU", nlon, nlev, nproma, ngptot, nblocks, plu);
  load_and_expand_2d(serializer, savepoint, "PLUDE", nlon, nlev, nproma, ngptot, nblocks, plude);
  load_and_expand_2d(serializer, savepoint, "PSNDE", nlon, nlev, nproma, ngptot, nblocks, psnde);
  load_and_expand_2d(serializer, savepoint, "PMFU", nlon, nlev, nproma, ngptot, nblocks, pmfu);
  load_and_expand_2d(serializer, savepoint, "PMFD", nlon, nlev, nproma, ngptot, nblocks, pmfd);
  load_and_expand_2d(serializer, savepoint, "PA", nlon, nlev, nproma, ngptot, nblocks, pa);
  load_and_expand_3d(serializer, savepoint, "PCLV", nlon, nlev, nclv, nproma, ngptot, nblocks, pclv);
  load_and_expand_2d(serializer, savepoint, "PSUPSAT", nlon, nlev, nproma, ngptot, nblocks, psupsat);


  *ptsphy = (float) serialboxMetainfoGetFloat64(metainfo, "PTSPHY");

  /* Populate global parameter values from meta-data */
  *rg = (float) serialboxMetainfoGetFloat64(metainfo, "RG");
  *rd = (float) serialboxMetainfoGetFloat64(metainfo, "RD");
  *rcpd = (float) serialboxMetainfoGetFloat64(metainfo, "RCPD");
  *retv = (float) serialboxMetainfoGetFloat64(metainfo, "RETV");
  *rlvtt = (float) serialboxMetainfoGetFloat64(metainfo, "RLVTT");
  *rlstt = (float) serialboxMetainfoGetFloat64(metainfo, "RLSTT");
  *rlmlt = (float) serialboxMetainfoGetFloat64(metainfo, "RLMLT");
  *rtt = (float) serialboxMetainfoGetFloat64(metainfo, "RTT");
  *rv = (float) serialboxMetainfoGetFloat64(metainfo, "RV");
  *r2es = (float) serialboxMetainfoGetFloat64(metainfo, "R2ES");
  *r3les = (float) serialboxMetainfoGetFloat64(metainfo, "R3LES");
  *r3ies = (float) serialboxMetainfoGetFloat64(metainfo, "R3IES");
  *r4les = (float) serialboxMetainfoGetFloat64(metainfo, "R4LES");
  *r4ies = (float) serialboxMetainfoGetFloat64(metainfo, "R4IES");
  *r5les = (float) serialboxMetainfoGetFloat64(metainfo, "R5LES");
  *r5ies = (float) serialboxMetainfoGetFloat64(metainfo, "R5IES");
  *r5alvcp = (float) serialboxMetainfoGetFloat64(metainfo, "R5ALVCP");
  *r5alscp = (float) serialboxMetainfoGetFloat64(metainfo, "R5ALSCP");
  *ralvdcp = (float) serialboxMetainfoGetFloat64(metainfo, "RALVDCP");
  *ralsdcp = (float) serialboxMetainfoGetFloat64(metainfo, "RALSDCP");
  *ralfdcp = (float) serialboxMetainfoGetFloat64(metainfo, "RALFDCP");
  *rtwat = (float) serialboxMetainfoGetFloat64(metainfo, "RTWAT");
  *rtice = (float) serialboxMetainfoGetFloat64(metainfo, "RTICE");
  *rticecu = (float) serialboxMetainfoGetFloat64(metainfo, "RTICECU");
  *rtwat_rtice_r = (float) serialboxMetainfoGetFloat64(metainfo, "RTWAT_RTICE_R");
  *rtwat_rticecu_r = (float) serialboxMetainfoGetFloat64(metainfo, "RTWAT_RTICECU_R");
  *rkoop1 = (float) serialboxMetainfoGetFloat64(metainfo, "RKOOP1");
  *rkoop2 = (float) serialboxMetainfoGetFloat64(metainfo, "RKOOP2");

  yrecldp->ramid = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RAMID");
  yrecldp->rcldiff = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLDIFF");
  yrecldp->rcldiff_convi = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLDIFF_CONVI");
  yrecldp->rclcrit = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLCRIT");
  yrecldp->rclcrit_sea = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLCRIT_SEA");
  yrecldp->rclcrit_land = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLCRIT_LAND");
  yrecldp->rkconv = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RKCONV");
  yrecldp->rprc1 = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RPRC1");
  yrecldp->rprc2 = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RPRC2");
  yrecldp->rcldmax = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLDMAX");
  yrecldp->rpecons = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RPECONS");
  yrecldp->rvrfactor = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RVRFACTOR");
  yrecldp->rprecrhmax = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RPRECRHMAX");
  yrecldp->rtaumel = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RTAUMEL");
  yrecldp->ramin = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RAMIN");
  yrecldp->rlmin = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RLMIN");
  yrecldp->rkooptau = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RKOOPTAU");

  yrecldp->rcldtopp = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLDTOPP");
  yrecldp->rlcritsnow = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RLCRITSNOW");
  yrecldp->rsnowlin1 = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RSNOWLIN1");
  yrecldp->rsnowlin2 = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RSNOWLIN2");
  yrecldp->ricehi1 = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RICEHI1");
  yrecldp->ricehi2 = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RICEHI2");
  yrecldp->riceinit = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RICEINIT");
  yrecldp->rvice = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RVICE");
  yrecldp->rvrain = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RVRAIN");
  yrecldp->rvsnow = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RVSNOW");
  yrecldp->rthomo = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RTHOMO");
  yrecldp->rcovpmin = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCOVPMIN");
  yrecldp->rccn = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCCN");
  yrecldp->rnice = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RNICE");
  yrecldp->rccnom = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCCNOM");
  yrecldp->rccnss = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCCNSS");
  yrecldp->rccnsu = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCCNSU");
  yrecldp->rcldtopcf = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLDTOPCF");
  yrecldp->rdepliqrefrate = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RDEPLIQREFRATE");
  yrecldp->rdepliqrefdepth = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RDEPLIQREFDEPTH");
  yrecldp->rcl_kkaac = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KKAac");
  yrecldp->rcl_kkbac = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KKBac");
  yrecldp->rcl_kkaau = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KKAau");
  yrecldp->rcl_kkbauq = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KKBauq");
  yrecldp->rcl_kkbaun = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KKBaun");
  yrecldp->rcl_kk_cloud_num_sea = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KK_cloud_num_sea");
  yrecldp->rcl_kk_cloud_num_land = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KK_cloud_num_land");
  yrecldp->rcl_ai = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_AI");
  yrecldp->rcl_bi = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_BI");
  yrecldp->rcl_ci = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CI");
  yrecldp->rcl_di = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_DI");
  yrecldp->rcl_x1i = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X1I");
  yrecldp->rcl_x2i = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X2I");
  yrecldp->rcl_x3i = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X3I");
  yrecldp->rcl_x4i = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X4I");
  yrecldp->rcl_const1i = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST1I");
  yrecldp->rcl_const2i = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST2I");
  yrecldp->rcl_const3i = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST3I");
  yrecldp->rcl_const4i = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST4I");
  yrecldp->rcl_const5i = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST5I");
  yrecldp->rcl_const6i = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST6I");
  yrecldp->rcl_apb1 = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_APB1");
  yrecldp->rcl_apb2 = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_APB2");
  yrecldp->rcl_apb3 = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_APB3");
  yrecldp->rcl_as = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_AS");
  yrecldp->rcl_bs = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_BS");
  yrecldp->rcl_cs = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CS");
  yrecldp->rcl_ds = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_DS");
  yrecldp->rcl_x1s = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X1S");
  yrecldp->rcl_x2s = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X2S");
  yrecldp->rcl_x3s = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X3S");
  yrecldp->rcl_x4s = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X4S");
  yrecldp->rcl_const1s = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST1S");
  yrecldp->rcl_const2s = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST2S");
  yrecldp->rcl_const3s = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST3S");
  yrecldp->rcl_const4s = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST4S");
  yrecldp->rcl_const5s = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST5S");
  yrecldp->rcl_const6s = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST6S");
  yrecldp->rcl_const7s = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST7S");
  yrecldp->rcl_const8s = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST8S");
  yrecldp->rdenswat = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RDENSWAT");
  yrecldp->rdensref = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RDENSREF");
  yrecldp->rcl_ar = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_AR");
  yrecldp->rcl_br = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_BR");
  yrecldp->rcl_cr = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CR");
  yrecldp->rcl_dr = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_DR");
  yrecldp->rcl_x1r = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X1R");
  yrecldp->rcl_x2r = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X2R");
  yrecldp->rcl_x4r = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X4R");
  yrecldp->rcl_ka273 = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KA273");
  yrecldp->rcl_cdenom1 = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CDENOM1");
  yrecldp->rcl_cdenom2 = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CDENOM2");
  yrecldp->rcl_cdenom3 = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CDENOM3");
  yrecldp->rcl_schmidt = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_SCHMIDT");
  yrecldp->rcl_dynvisc = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_DYNVISC");
  yrecldp->rcl_const1r = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST1R");
  yrecldp->rcl_const2r = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST2R");
  yrecldp->rcl_const3r = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST3R");
  yrecldp->rcl_const4r = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST4R");
  yrecldp->rcl_fac1 = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_FAC1");
  yrecldp->rcl_fac2 = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_FAC2");
  yrecldp->rcl_const5r = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST5R");
  yrecldp->rcl_const6r = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST6R");
  yrecldp->rcl_fzrab = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_FZRAB");
  yrecldp->rcl_fzrbb = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_FZRBB");
  yrecldp->lcldextra = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LCLDEXTRA");
  yrecldp->lcldbudget = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LCLDBUDGET");
  yrecldp->nssopt = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NSSOPT");
  yrecldp->ncldtop = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NCLDTOP");
  yrecldp->naeclbc = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAECLBC");
  yrecldp->naecldu = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAECLDU");
  yrecldp->naeclom = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAECLOM");
  yrecldp->naeclss = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAECLSS");
  yrecldp->naeclsu = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAECLSU");
  yrecldp->nclddiag = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NCLDDIAG");
  yrecldp->naercld = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAERCLD");
  yrecldp->laerliqautolsp = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERLIQAUTOLSP");
  yrecldp->laerliqautocp = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERLIQAUTOCP");
  yrecldp->laerliqautocpb = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERLIQAUTOCPB");
  yrecldp->laerliqcoll = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERLIQCOLL");
  yrecldp->laericesed = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERICESED");
  yrecldp->laericeauto = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERICEAUTO");
  yrecldp->nshapep = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NSHAPEP");
  yrecldp->nshapeq = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NSHAPEQ");
  yrecldp->nbeta = (float) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NBETA");

  serialboxSerializerDestroySavepointVector(savepoints, 1);
  serialboxSerializerDestroy(serializer);
}

/* Read reference result into memory */
void load_reference(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma,
		    float *plude, float *pcovptot, float *prainfrac_toprfz, float *pfsqlf, float *pfsqif,
		    float *pfcqlng, float *pfcqnng, float *pfsqrf, float *pfsqsf, float *pfcqrng, float *pfcqsng,
		    float *pfsqltur, float *pfsqitur, float *pfplsl, float *pfplsn, float *pfhpsl, float *pfhpsn,
		    float *tend_loc_a, float *tend_loc_q, float *tend_loc_t, float *tend_loc_cld)
{
  serialboxSerializer_t* serializer = serialboxSerializerCreate(Read, "./data", "reference", "Binary");
  serialboxMetainfo_t* metainfo = serialboxSerializerGetGlobalMetainfo(serializer);
  serialboxSavepoint_t** savepoints = serialboxSerializerGetSavepointVector(serializer);
  serialboxSavepoint_t* savepoint = savepoints[0];

  int nblocks = (ngptot / nproma) + min(ngptot % nproma, 1);

  load_and_expand_2d(serializer, savepoint, "PLUDE", nlon, nlev, nproma, ngptot, nblocks, plude);
  load_and_expand_2d(serializer, savepoint, "PCOVPTOT", nlon, nlev, nproma, ngptot, nblocks, pcovptot);
  load_and_expand_1d(serializer, savepoint, "PRAINFRAC_TOPRFZ", nlon, nproma, ngptot, nblocks, prainfrac_toprfz);
  load_and_expand_2d(serializer, savepoint, "PFSQLF", nlon, nlev+1, nproma, ngptot, nblocks, pfsqlf);
  load_and_expand_2d(serializer, savepoint, "PFSQIF", nlon, nlev+1, nproma, ngptot, nblocks, pfsqif);
  load_and_expand_2d(serializer, savepoint, "PFCQLNG", nlon, nlev+1, nproma, ngptot, nblocks, pfcqlng);
  load_and_expand_2d(serializer, savepoint, "PFCQNNG", nlon, nlev+1, nproma, ngptot, nblocks, pfcqnng);
  load_and_expand_2d(serializer, savepoint, "PFSQRF", nlon, nlev+1, nproma, ngptot, nblocks, pfsqrf);
  load_and_expand_2d(serializer, savepoint, "PFSQSF", nlon, nlev+1, nproma, ngptot, nblocks, pfsqsf);
  load_and_expand_2d(serializer, savepoint, "PFCQRNG", nlon, nlev+1, nproma, ngptot, nblocks, pfcqrng);
  load_and_expand_2d(serializer, savepoint, "PFCQSNG", nlon, nlev+1, nproma, ngptot, nblocks, pfcqsng);
  load_and_expand_2d(serializer, savepoint, "PFSQLTUR", nlon, nlev+1, nproma, ngptot, nblocks, pfsqltur);
  load_and_expand_2d(serializer, savepoint, "PFSQITUR", nlon, nlev+1, nproma, ngptot, nblocks, pfsqitur);
  load_and_expand_2d(serializer, savepoint, "PFPLSL", nlon, nlev+1, nproma, ngptot, nblocks, pfplsl);
  load_and_expand_2d(serializer, savepoint, "PFPLSN", nlon, nlev+1, nproma, ngptot, nblocks, pfplsn);
  load_and_expand_2d(serializer, savepoint, "PFHPSL", nlon, nlev+1, nproma, ngptot, nblocks, pfhpsl);
  load_and_expand_2d(serializer, savepoint, "PFHPSN", nlon, nlev+1, nproma, ngptot, nblocks, pfhpsn);
  load_and_expand_2d(serializer, savepoint, "TENDENCY_LOC_T", nlon, nlev, nproma, ngptot, nblocks, tend_loc_t);
  load_and_expand_2d(serializer, savepoint, "TENDENCY_LOC_Q", nlon, nlev, nproma, ngptot, nblocks, tend_loc_q);
  load_and_expand_2d(serializer, savepoint, "TENDENCY_LOC_A", nlon, nlev, nproma, ngptot, nblocks, tend_loc_a);
  load_and_expand_3d(serializer, savepoint, "TENDENCY_LOC_CLD", nlon, nlev, nclv, nproma, ngptot, nblocks, tend_loc_cld);

  serialboxSerializerDestroySavepointVector(savepoints, 1);
  serialboxSerializerDestroy(serializer);
}
