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
#include <iostream>

#include <math.h>
#ifdef HAVE_SERIALBOX
#include "serialbox-c/Serialbox.h"
#endif
#ifdef HAVE_HDF5
#include "hdf5.h"
#define INPUT_FILE "input.h5"
#define REFERENCE_FILE "reference.h5"
#endif

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))

#ifdef HAVE_HDF5
void read_hdf5_int(hid_t file_id, const char *name, int *field) {
  hid_t dataset_id;
  herr_t  status;
  dataset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, field);
  status = H5Dclose(dataset_id);
}

void read_hdf5(hid_t file_id, const char *name, dtype *field) {
  hid_t dataset_id;
  herr_t  status;
  double dbl_field;
  dataset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, &dbl_field);
  *field = (dtype)dbl_field;
  status = H5Dclose(dataset_id);
}
#endif

/* Query sizes and dimensions of state arrays */
void query_state(int *klon, int *klev)
{
#ifdef HAVE_SERIALBOX
  serialboxSerializer_t* serializer = serialboxSerializerCreate(Read, "./data", "input", "Binary");
  serialboxMetainfo_t* globalMetainfo = serialboxSerializerGetGlobalMetainfo(serializer);

  *klon = serialboxMetainfoGetInt32(globalMetainfo, "KLON");
  *klev = serialboxMetainfoGetInt32(globalMetainfo, "KLEV");

  serialboxMetainfoDestroy(globalMetainfo);
  serialboxSerializerDestroy(serializer);
#endif
#ifdef HAVE_HDF5
  hid_t file_id, dataset_id;
  herr_t  status;
  file_id = H5Fopen(INPUT_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

  read_hdf5_int(file_id, "/KLEV", klev);
  read_hdf5_int(file_id, "/KLON", klon);

  status = H5Fclose(file_id);
#endif
}

void expand_1d(dtype *buffer, dtype *field_in, int nlon, int nproma, int ngptot, int nblocks)
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


void expand_2d(dtype *buffer_in, dtype *field_in, int nlon, int nlev, int nproma, int ngptot, int nblocks)
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

void expand_3d(dtype *buffer_in, dtype *field_in, int nlon, int nlev, int nclv, int nproma, int ngptot, int nblocks)
{
  int b, l, c, i, buf_start_idx, buf_idx;

  #pragma omp parallel for default(shared) private(b, buf_start_idx, buf_idx, l, i, c)
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

#ifdef HAVE_SERIALBOX
void load_and_expand_1d(serialboxSerializer_t *serializer, serialboxSavepoint_t* savepoint,
    const char *name, int nlon, int nproma, int ngptot, int nblocks, dtype *field)
{
  double dbl_buffer[nlon];
  dtype buffer[nlon];
  int strides[1] = {1};

  serialboxSerializerRead(serializer, name, savepoint, dbl_buffer, strides, 1);
  for (int i=0;i<nlon;i++) {
    buffer[i] = (dtype) dbl_buffer[i];
  }
  expand_1d((dtype *)buffer, field, nlon, nproma, ngptot, nblocks);
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
    const char *name, int nlon, int nlev, int nproma, int ngptot, int nblocks, dtype *field)
{
  double dbl_buffer[nlev][nlon];
  dtype buffer[nlev][nlon];
  int strides[2] = {1, nlon};

  serialboxSerializerRead(serializer, name, savepoint, dbl_buffer, strides, 2);
  for (int j=0;j<nlev;j++) {
  for (int i=0;i<nlon;i++) {
    buffer[j][i] = (dtype) dbl_buffer[j][i];
  }
  }
  expand_2d((dtype *)buffer, field, nlon, nlev, nproma, ngptot, nblocks);
}

void load_and_expand_3d(serialboxSerializer_t *serializer, serialboxSavepoint_t* savepoint,
    const char *name, int nlon, int nlev, int nclv, int nproma, int ngptot, int nblocks, dtype *field)
{
  double dbl_buffer[nclv][nlev][nlon];
  dtype buffer[nclv][nlev][nlon];
  int strides[3] = {1, nlon, nlev*nlon};

  serialboxSerializerRead(serializer, name, savepoint, dbl_buffer, strides, 3);
  for (int k=0;k<nclv;k++) {
  for (int j=0;j<nlev;j++) {
  for (int i=0;i<nlon;i++) {
    buffer[k][j][i] = (dtype) dbl_buffer[k][j][i];
  }
  }
  }
  expand_3d((dtype *)buffer, field, nlon, nlev, nclv, nproma, ngptot, nblocks);
}
#endif

#if HAVE_HDF5
void load_and_expand_1d(hid_t file_id, const char *name, int nlon, int nproma, int ngptot, int nblocks, dtype *field)
{
  dtype buffer[nlon];
  double dbl_buffer[nlon];
  int strides[1] = {1};
  hid_t dataset_id;
  dataset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
  herr_t  status;
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dbl_buffer);
  status = H5Dclose(dataset_id);
  for (int i=0;i<nlon;i++) {
    buffer[i] = (dtype) dbl_buffer[i];
  }
  expand_1d((dtype *)buffer, field, nlon, nproma, ngptot, nblocks);
}

void load_and_expand_1d_int(hid_t file_id, const char *name, int nlon, int nproma, int ngptot, int nblocks, int *field)
{
  int buffer[nlon];
  int strides[1] = {1};
  hid_t dataset_id;
  dataset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
  herr_t  status;
  status = H5Dread(dataset_id, H5T_NATIVE_INT, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  status = H5Dclose(dataset_id);
  expand_1d_int((int *)buffer, field, nlon, nproma, ngptot, nblocks);
}

void load_and_expand_2d(hid_t file_id, const char *name, int nlon, int nlev, int nproma, int ngptot, int nblocks, dtype *field)
{
  dtype buffer[nlev][nlon];
  double dbl_buffer[nlev][nlon];
  int strides[2] = {1, nlon};
  hid_t dataset_id;
  dataset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
  herr_t  status;
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dbl_buffer);
  status = H5Dclose(dataset_id);
  for (int j=0;j<nlev;j++) {
  for (int i=0;i<nlon;i++) {
    buffer[j][i] = (dtype) dbl_buffer[j][i];
  }
  }
  expand_2d((dtype *)buffer, field, nlon, nlev, nproma, ngptot, nblocks);
}

void load_and_expand_3d(hid_t file_id, const char *name, int nlon, int nlev, int nclv, int nproma, int ngptot, int nblocks, dtype *field)
{
  dtype buffer[nclv][nlev][nlon];
  double dbl_buffer[nclv][nlev][nlon];
  int strides[3] = {1, nlon, nlev*nlon};
  hid_t dataset_id;
  dataset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
  herr_t  status;
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, dbl_buffer);
  status = H5Dclose(dataset_id);
  for (int k=0;k<nclv;k++) {
  for (int j=0;j<nlev;j++) {
  for (int i=0;i<nlon;i++) {
    buffer[k][j][i] = (dtype) dbl_buffer[k][j][i];
  }
  }
  }
  expand_3d((dtype *)buffer, field, nlon, nlev, nclv, nproma, ngptot, nblocks);
}
#endif


/* Read input state into memory */
void load_state(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma,
    dtype* ptsphy, dtype* plcrit_aer, dtype* picrit_aer,
    dtype* pre_ice, dtype* pccn, dtype* pnice, dtype* pt, dtype* pq,
    dtype* tend_cml_t, dtype* tend_cml_q, dtype* tend_cml_a, dtype* tend_cml_cld,
    dtype* tend_tmp_t, dtype* tend_tmp_q, dtype* tend_tmp_a, dtype* tend_tmp_cld,
    dtype* pvfa, dtype* pvfl, dtype* pvfi, dtype* pdyna, dtype* pdynl, dtype* pdyni,
    dtype* phrsw, dtype* phrlw, dtype* pvervel, dtype* pap, dtype* paph, dtype* plsm,
    int* ktype, dtype* plu, dtype* plude, dtype* psnde, dtype* pmfu,
    dtype* pmfd, dtype* pa, dtype* pclv, dtype* psupsat, struct TECLDP* yrecldp,
    dtype* rg, dtype* rd, dtype* rcpd, dtype* retv, dtype* rlvtt, dtype* rlstt,
    dtype* rlmlt, dtype* rtt, dtype* rv, dtype* r2es, dtype* r3les, dtype* r3ies,
    dtype* r4les, dtype* r4ies, dtype* r5les, dtype* r5ies, dtype* r5alvcp, dtype* r5alscp,
    dtype* ralvdcp, dtype* ralsdcp, dtype* ralfdcp, dtype* rtwat,
    dtype* rtice, dtype* rticecu, dtype* rtwat_rtice_r, dtype *rtwat_rticecu_r,
    dtype* rkoop1, dtype* rkoop2 )
{

  int nblocks = (ngptot / nproma) + min(ngptot % nproma, 1);

#ifdef HAVE_SERIALBOX
  serialboxSerializer_t* serializer = serialboxSerializerCreate(Read, "./data", "input", "Binary");
  serialboxMetainfo_t* metainfo = serialboxSerializerGetGlobalMetainfo(serializer);
  serialboxSavepoint_t** savepoints = serialboxSerializerGetSavepointVector(serializer);
  serialboxSavepoint_t* savepoint = savepoints[0];

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


  *ptsphy = (dtype) serialboxMetainfoGetFloat64(metainfo, "PTSPHY");

  /* Populate global parameter values from meta-data */
  *rg = (dtype) serialboxMetainfoGetFloat64(metainfo, "RG");
  *rd = (dtype) serialboxMetainfoGetFloat64(metainfo, "RD");
  *rcpd = (dtype) serialboxMetainfoGetFloat64(metainfo, "RCPD");
  *retv = (dtype) serialboxMetainfoGetFloat64(metainfo, "RETV");
  *rlvtt = (dtype) serialboxMetainfoGetFloat64(metainfo, "RLVTT");
  *rlstt = (dtype) serialboxMetainfoGetFloat64(metainfo, "RLSTT");
  *rlmlt = (dtype) serialboxMetainfoGetFloat64(metainfo, "RLMLT");
  *rtt = (dtype) serialboxMetainfoGetFloat64(metainfo, "RTT");
  *rv = (dtype) serialboxMetainfoGetFloat64(metainfo, "RV");
  *r2es = (dtype) serialboxMetainfoGetFloat64(metainfo, "R2ES");
  *r3les = (dtype) serialboxMetainfoGetFloat64(metainfo, "R3LES");
  *r3ies = (dtype) serialboxMetainfoGetFloat64(metainfo, "R3IES");
  *r4les = (dtype) serialboxMetainfoGetFloat64(metainfo, "R4LES");
  *r4ies = (dtype) serialboxMetainfoGetFloat64(metainfo, "R4IES");
  *r5les = (dtype) serialboxMetainfoGetFloat64(metainfo, "R5LES");
  *r5ies = (dtype) serialboxMetainfoGetFloat64(metainfo, "R5IES");
  *r5alvcp = (dtype) serialboxMetainfoGetFloat64(metainfo, "R5ALVCP");
  *r5alscp = (dtype) serialboxMetainfoGetFloat64(metainfo, "R5ALSCP");
  *ralvdcp = (dtype) serialboxMetainfoGetFloat64(metainfo, "RALVDCP");
  *ralsdcp = (dtype) serialboxMetainfoGetFloat64(metainfo, "RALSDCP");
  *ralfdcp = (dtype) serialboxMetainfoGetFloat64(metainfo, "RALFDCP");
  *rtwat = (dtype) serialboxMetainfoGetFloat64(metainfo, "RTWAT");
  *rtice = (dtype) serialboxMetainfoGetFloat64(metainfo, "RTICE");
  *rticecu = (dtype) serialboxMetainfoGetFloat64(metainfo, "RTICECU");
  *rtwat_rtice_r = (dtype) serialboxMetainfoGetFloat64(metainfo, "RTWAT_RTICE_R");
  *rtwat_rticecu_r = (dtype) serialboxMetainfoGetFloat64(metainfo, "RTWAT_RTICECU_R");
  *rkoop1 = (dtype) serialboxMetainfoGetFloat64(metainfo, "RKOOP1");
  *rkoop2 = (dtype) serialboxMetainfoGetFloat64(metainfo, "RKOOP2");

  yrecldp->ramid = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RAMID");
  yrecldp->rcldiff = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLDIFF");
  yrecldp->rcldiff_convi = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLDIFF_CONVI");
  yrecldp->rclcrit = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLCRIT");
  yrecldp->rclcrit_sea = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLCRIT_SEA");
  yrecldp->rclcrit_land = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLCRIT_LAND");
  yrecldp->rkconv = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RKCONV");
  yrecldp->rprc1 = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RPRC1");
  yrecldp->rprc2 = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RPRC2");
  yrecldp->rcldmax = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLDMAX");
  yrecldp->rpecons = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RPECONS");
  yrecldp->rvrfactor = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RVRFACTOR");
  yrecldp->rprecrhmax = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RPRECRHMAX");
  yrecldp->rtaumel = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RTAUMEL");
  yrecldp->ramin = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RAMIN");
  yrecldp->rlmin = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RLMIN");
  yrecldp->rkooptau = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RKOOPTAU");

  yrecldp->rcldtopp = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLDTOPP");
  yrecldp->rlcritsnow = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RLCRITSNOW");
  yrecldp->rsnowlin1 = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RSNOWLIN1");
  yrecldp->rsnowlin2 = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RSNOWLIN2");
  yrecldp->ricehi1 = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RICEHI1");
  yrecldp->ricehi2 = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RICEHI2");
  yrecldp->riceinit = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RICEINIT");
  yrecldp->rvice = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RVICE");
  yrecldp->rvrain = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RVRAIN");
  yrecldp->rvsnow = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RVSNOW");
  yrecldp->rthomo = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RTHOMO");
  yrecldp->rcovpmin = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCOVPMIN");
  yrecldp->rccn = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCCN");
  yrecldp->rnice = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RNICE");
  yrecldp->rccnom = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCCNOM");
  yrecldp->rccnss = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCCNSS");
  yrecldp->rccnsu = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCCNSU");
  yrecldp->rcldtopcf = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLDTOPCF");
  yrecldp->rdepliqrefrate = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RDEPLIQREFRATE");
  yrecldp->rdepliqrefdepth = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RDEPLIQREFDEPTH");
  yrecldp->rcl_kkaac = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KKAac");
  yrecldp->rcl_kkbac = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KKBac");
  yrecldp->rcl_kkaau = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KKAau");
  yrecldp->rcl_kkbauq = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KKBauq");
  yrecldp->rcl_kkbaun = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KKBaun");
  yrecldp->rcl_kk_cloud_num_sea = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KK_cloud_num_sea");
  yrecldp->rcl_kk_cloud_num_land = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KK_cloud_num_land");
  yrecldp->rcl_ai = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_AI");
  yrecldp->rcl_bi = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_BI");
  yrecldp->rcl_ci = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CI");
  yrecldp->rcl_di = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_DI");
  yrecldp->rcl_x1i = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X1I");
  yrecldp->rcl_x2i = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X2I");
  yrecldp->rcl_x3i = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X3I");
  yrecldp->rcl_x4i = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X4I");
  yrecldp->rcl_const1i = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST1I");
  yrecldp->rcl_const2i = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST2I");
  yrecldp->rcl_const3i = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST3I");
  yrecldp->rcl_const4i = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST4I");
  yrecldp->rcl_const5i = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST5I");
  yrecldp->rcl_const6i = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST6I");
  yrecldp->rcl_apb1 = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_APB1");
  yrecldp->rcl_apb2 = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_APB2");
  yrecldp->rcl_apb3 = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_APB3");
  yrecldp->rcl_as = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_AS");
  yrecldp->rcl_bs = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_BS");
  yrecldp->rcl_cs = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CS");
  yrecldp->rcl_ds = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_DS");
  yrecldp->rcl_x1s = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X1S");
  yrecldp->rcl_x2s = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X2S");
  yrecldp->rcl_x3s = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X3S");
  yrecldp->rcl_x4s = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X4S");
  yrecldp->rcl_const1s = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST1S");
  yrecldp->rcl_const2s = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST2S");
  yrecldp->rcl_const3s = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST3S");
  yrecldp->rcl_const4s = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST4S");
  yrecldp->rcl_const5s = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST5S");
  yrecldp->rcl_const6s = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST6S");
  yrecldp->rcl_const7s = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST7S");
  yrecldp->rcl_const8s = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST8S");
  yrecldp->rdenswat = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RDENSWAT");
  yrecldp->rdensref = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RDENSREF");
  yrecldp->rcl_ar = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_AR");
  yrecldp->rcl_br = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_BR");
  yrecldp->rcl_cr = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CR");
  yrecldp->rcl_dr = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_DR");
  yrecldp->rcl_x1r = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X1R");
  yrecldp->rcl_x2r = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X2R");
  yrecldp->rcl_x4r = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X4R");
  yrecldp->rcl_ka273 = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KA273");
  yrecldp->rcl_cdenom1 = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CDENOM1");
  yrecldp->rcl_cdenom2 = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CDENOM2");
  yrecldp->rcl_cdenom3 = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CDENOM3");
  yrecldp->rcl_schmidt = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_SCHMIDT");
  yrecldp->rcl_dynvisc = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_DYNVISC");
  yrecldp->rcl_const1r = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST1R");
  yrecldp->rcl_const2r = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST2R");
  yrecldp->rcl_const3r = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST3R");
  yrecldp->rcl_const4r = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST4R");
  yrecldp->rcl_fac1 = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_FAC1");
  yrecldp->rcl_fac2 = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_FAC2");
  yrecldp->rcl_const5r = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST5R");
  yrecldp->rcl_const6r = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST6R");
  yrecldp->rcl_fzrab = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_FZRAB");
  yrecldp->rcl_fzrbb = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_FZRBB");
  yrecldp->lcldextra = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LCLDEXTRA");
  yrecldp->lcldbudget = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LCLDBUDGET");
  yrecldp->nssopt = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NSSOPT");
  yrecldp->ncldtop = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NCLDTOP");
  yrecldp->naeclbc = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAECLBC");
  yrecldp->naecldu = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAECLDU");
  yrecldp->naeclom = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAECLOM");
  yrecldp->naeclss = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAECLSS");
  yrecldp->naeclsu = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAECLSU");
  yrecldp->nclddiag = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NCLDDIAG");
  yrecldp->naercld = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAERCLD");
  yrecldp->laerliqautolsp = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERLIQAUTOLSP");
  yrecldp->laerliqautocp = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERLIQAUTOCP");
  yrecldp->laerliqautocpb = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERLIQAUTOCPB");
  yrecldp->laerliqcoll = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERLIQCOLL");
  yrecldp->laericesed = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERICESED");
  yrecldp->laericeauto = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERICEAUTO");
  yrecldp->nshapep = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NSHAPEP");
  yrecldp->nshapeq = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NSHAPEQ");
  yrecldp->nbeta = (dtype) serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NBETA");

  serialboxSerializerDestroySavepointVector(savepoints, 1);
  serialboxSerializerDestroy(serializer);
#endif

#ifdef HAVE_HDF5

  hid_t file_id, dataset_id;
  herr_t  status;
  file_id = H5Fopen(INPUT_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

  load_and_expand_2d(file_id, "PLCRIT_AER", nlon, nlev, nproma, ngptot, nblocks, plcrit_aer);
  load_and_expand_2d(file_id, "PICRIT_AER", nlon, nlev, nproma, ngptot, nblocks, picrit_aer);
  load_and_expand_2d(file_id, "PRE_ICE", nlon, nlev, nproma, ngptot, nblocks, pre_ice);
  load_and_expand_2d(file_id, "PCCN", nlon, nlev, nproma, ngptot, nblocks, pccn);
  load_and_expand_2d(file_id, "PNICE", nlon, nlev, nproma, ngptot, nblocks, pnice);
  load_and_expand_2d(file_id, "PT", nlon, nlev, nproma, ngptot, nblocks, pt);
  load_and_expand_2d(file_id, "PQ", nlon, nlev, nproma, ngptot, nblocks, pq);
  load_and_expand_2d(file_id, "TENDENCY_CML_T", nlon, nlev, nproma, ngptot, nblocks, tend_cml_t);
  load_and_expand_2d(file_id, "TENDENCY_CML_Q", nlon, nlev, nproma, ngptot, nblocks, tend_cml_q);
  load_and_expand_2d(file_id, "TENDENCY_CML_A", nlon, nlev, nproma, ngptot, nblocks, tend_cml_a);
  load_and_expand_3d(file_id, "TENDENCY_CML_CLD", nlon, nlev, nclv, nproma, ngptot, nblocks, tend_cml_cld);
  load_and_expand_2d(file_id, "TENDENCY_TMP_T", nlon, nlev, nproma, ngptot, nblocks, tend_tmp_t);
  load_and_expand_2d(file_id, "TENDENCY_TMP_Q", nlon, nlev, nproma, ngptot, nblocks, tend_tmp_q);
  load_and_expand_2d(file_id, "TENDENCY_TMP_A", nlon, nlev, nproma, ngptot, nblocks, tend_tmp_a);
  load_and_expand_3d(file_id, "TENDENCY_TMP_CLD", nlon, nlev, nclv, nproma, ngptot, nblocks, tend_tmp_cld);
  load_and_expand_2d(file_id, "PVFA", nlon, nlev, nproma, ngptot, nblocks, pvfa);
  load_and_expand_2d(file_id, "PVFL", nlon, nlev, nproma, ngptot, nblocks, pvfl);
  load_and_expand_2d(file_id, "PVFI", nlon, nlev, nproma, ngptot, nblocks, pvfi);
  load_and_expand_2d(file_id, "PDYNA", nlon, nlev, nproma, ngptot, nblocks, pdyna);
  load_and_expand_2d(file_id, "PDYNL", nlon, nlev, nproma, ngptot, nblocks, pdynl);
  load_and_expand_2d(file_id, "PDYNI", nlon, nlev, nproma, ngptot, nblocks, pdyni);
  load_and_expand_2d(file_id, "PHRSW", nlon, nlev, nproma, ngptot, nblocks, phrsw);
  load_and_expand_2d(file_id, "PHRLW", nlon, nlev, nproma, ngptot, nblocks, phrlw);
  load_and_expand_2d(file_id, "PVERVEL", nlon, nlev, nproma, ngptot, nblocks, pvervel);
  load_and_expand_2d(file_id, "PAP", nlon, nlev, nproma, ngptot, nblocks, pap);
  load_and_expand_2d(file_id, "PAPH", nlon, nlev+1, nproma, ngptot, nblocks, paph);
  load_and_expand_1d(file_id, "PLSM", nlon, nproma, ngptot, nblocks, plsm);
  load_and_expand_1d_int(file_id, "KTYPE", nlon, nproma, ngptot, nblocks, ktype);
  load_and_expand_2d(file_id, "PLU", nlon, nlev, nproma, ngptot, nblocks, plu);
  load_and_expand_2d(file_id, "PLUDE", nlon, nlev, nproma, ngptot, nblocks, plude);
  load_and_expand_2d(file_id, "PSNDE", nlon, nlev, nproma, ngptot, nblocks, psnde);
  load_and_expand_2d(file_id, "PMFU", nlon, nlev, nproma, ngptot, nblocks, pmfu);
  load_and_expand_2d(file_id, "PMFD", nlon, nlev, nproma, ngptot, nblocks, pmfd);
  load_and_expand_2d(file_id, "PA", nlon, nlev, nproma, ngptot, nblocks, pa);
  load_and_expand_3d(file_id, "PCLV", nlon, nlev, nclv, nproma, ngptot, nblocks, pclv);
  load_and_expand_2d(file_id, "PSUPSAT", nlon, nlev, nproma, ngptot, nblocks, psupsat);

  read_hdf5(file_id, "/PTSPHY", ptsphy);

  read_hdf5(file_id, "/RG", rg);
  read_hdf5(file_id, "/RD", rd);
  read_hdf5(file_id, "/RCPD", rcpd);
  read_hdf5(file_id, "/RETV", retv);
  read_hdf5(file_id, "/RLVTT", rlvtt);
  read_hdf5(file_id, "/RLSTT", rlstt);
  read_hdf5(file_id, "/RLMLT", rlmlt);
  read_hdf5(file_id, "/RTT", rtt);
  read_hdf5(file_id, "/RV", rv);
  read_hdf5(file_id, "/R2ES", r2es);
  read_hdf5(file_id, "/R3LES", r3les);
  read_hdf5(file_id, "/R3IES", r3ies);
  read_hdf5(file_id, "/R4LES", r4les);
  read_hdf5(file_id, "/R4IES", r4ies);
  read_hdf5(file_id, "/R5LES", r5les);
  read_hdf5(file_id, "/R5IES", r5ies);
  read_hdf5(file_id, "/R5ALVCP", r5alvcp);
  read_hdf5(file_id, "/R5ALSCP", r5alscp);
  read_hdf5(file_id, "/RALVDCP", ralvdcp);
  read_hdf5(file_id, "/RALSDCP", ralsdcp);
  read_hdf5(file_id, "/RALFDCP", ralfdcp);
  read_hdf5(file_id, "/RTWAT", rtwat);
  read_hdf5(file_id, "/RTICE", rtice);
  read_hdf5(file_id, "/RTICECU", rticecu);
  read_hdf5(file_id, "/RTWAT_RTICE_R", rtwat_rtice_r);
  read_hdf5(file_id, "/RTWAT_RTICECU_R", rtwat_rticecu_r);
  read_hdf5(file_id, "/RKOOP1", rkoop1);
  read_hdf5(file_id, "/RKOOP2", rkoop2);

  read_hdf5(file_id, "/YRECLDP_RAMID", &yrecldp->ramid);
  read_hdf5(file_id, "/YRECLDP_RCLDIFF", &yrecldp->rcldiff);
  read_hdf5(file_id, "/YRECLDP_RCLDIFF_CONVI", &yrecldp->rcldiff_convi);
  read_hdf5(file_id, "/YRECLDP_RCLCRIT", &yrecldp->rclcrit);
  read_hdf5(file_id, "/YRECLDP_RCLCRIT_SEA", &yrecldp->rclcrit_sea);
  read_hdf5(file_id, "/YRECLDP_RCLCRIT_LAND", &yrecldp->rclcrit_land);
  read_hdf5(file_id, "/YRECLDP_RKCONV", &yrecldp->rkconv);
  read_hdf5(file_id, "/YRECLDP_RPRC1", &yrecldp->rprc1);
  read_hdf5(file_id, "/YRECLDP_RPRC2", &yrecldp->rprc2);
  read_hdf5(file_id, "/YRECLDP_RCLDMAX", &yrecldp->rcldmax);
  read_hdf5(file_id, "/YRECLDP_RPECONS", &yrecldp->rpecons);
  read_hdf5(file_id, "/YRECLDP_RVRFACTOR", &yrecldp->rvrfactor);
  read_hdf5(file_id, "/YRECLDP_RPRECRHMAX", &yrecldp->rprecrhmax);
  read_hdf5(file_id, "/YRECLDP_RTAUMEL", &yrecldp->rtaumel);
  read_hdf5(file_id, "/YRECLDP_RAMIN", &yrecldp->ramin);
  read_hdf5(file_id, "/YRECLDP_RLMIN", &yrecldp->rlmin);
  read_hdf5(file_id, "/YRECLDP_RKOOPTAU", &yrecldp->rkooptau);
  read_hdf5(file_id, "/YRECLDP_RCLDTOPP", &yrecldp->rcldtopp);
  read_hdf5(file_id, "/YRECLDP_RLCRITSNOW", &yrecldp->rlcritsnow);
  read_hdf5(file_id, "/YRECLDP_RSNOWLIN1", &yrecldp->rsnowlin1);
  read_hdf5(file_id, "/YRECLDP_RSNOWLIN2", &yrecldp->rsnowlin2);
  read_hdf5(file_id, "/YRECLDP_RICEHI1", &yrecldp->ricehi1);
  read_hdf5(file_id, "/YRECLDP_RICEHI2", &yrecldp->ricehi2);
  read_hdf5(file_id, "/YRECLDP_RICEINIT", &yrecldp->riceinit);
  read_hdf5(file_id, "/YRECLDP_RVICE", &yrecldp->rvice);
  read_hdf5(file_id, "/YRECLDP_RVRAIN", &yrecldp->rvrain);
  read_hdf5(file_id, "/YRECLDP_RVSNOW", &yrecldp->rvsnow);
  read_hdf5(file_id, "/YRECLDP_RTHOMO", &yrecldp->rthomo);
  read_hdf5(file_id, "/YRECLDP_RCOVPMIN", &yrecldp->rcovpmin);
  read_hdf5(file_id, "/YRECLDP_RCCN", &yrecldp->rccn);
  read_hdf5(file_id, "/YRECLDP_RNICE", &yrecldp->rnice);
  read_hdf5(file_id, "/YRECLDP_RCCNOM", &yrecldp->rccnom);
  read_hdf5(file_id, "/YRECLDP_RCCNSS", &yrecldp->rccnss);
  read_hdf5(file_id, "/YRECLDP_RCCNSU", &yrecldp->rccnsu);
  read_hdf5(file_id, "/YRECLDP_RCLDTOPCF", &yrecldp->rcldtopcf);
  read_hdf5(file_id, "/YRECLDP_RDEPLIQREFRATE", &yrecldp->rdepliqrefrate);
  read_hdf5(file_id, "/YRECLDP_RDEPLIQREFDEPTH", &yrecldp->rdepliqrefdepth);
  read_hdf5(file_id, "/YRECLDP_RCL_KKAac", &yrecldp->rcl_kkaac);
  read_hdf5(file_id, "/YRECLDP_RCL_KKBac", &yrecldp->rcl_kkbac);
  read_hdf5(file_id, "/YRECLDP_RCL_KKAau", &yrecldp->rcl_kkaau);
  read_hdf5(file_id, "/YRECLDP_RCL_KKBauq", &yrecldp->rcl_kkbauq);
  read_hdf5(file_id, "/YRECLDP_RCL_KKBaun", &yrecldp->rcl_kkbaun);
  read_hdf5(file_id, "/YRECLDP_RCL_KK_cloud_num_sea", &yrecldp->rcl_kk_cloud_num_sea);
  read_hdf5(file_id, "/YRECLDP_RCL_KK_cloud_num_land", &yrecldp->rcl_kk_cloud_num_land);
  read_hdf5(file_id, "/YRECLDP_RCL_AI", &yrecldp->rcl_ai);
  read_hdf5(file_id, "/YRECLDP_RCL_BI", &yrecldp->rcl_bi);
  read_hdf5(file_id, "/YRECLDP_RCL_CI", &yrecldp->rcl_ci);
  read_hdf5(file_id, "/YRECLDP_RCL_DI", &yrecldp->rcl_di);
  read_hdf5(file_id, "/YRECLDP_RCL_X1I", &yrecldp->rcl_x1i);
  read_hdf5(file_id, "/YRECLDP_RCL_X2I", &yrecldp->rcl_x2i);
  read_hdf5(file_id, "/YRECLDP_RCL_X3I", &yrecldp->rcl_x3i);
  read_hdf5(file_id, "/YRECLDP_RCL_X4I", &yrecldp->rcl_x4i);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST1I", &yrecldp->rcl_const1i);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST2I", &yrecldp->rcl_const2i);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST3I", &yrecldp->rcl_const3i);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST4I", &yrecldp->rcl_const4i);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST5I", &yrecldp->rcl_const5i);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST6I", &yrecldp->rcl_const6i);
  read_hdf5(file_id, "/YRECLDP_RCL_APB1", &yrecldp->rcl_apb1);
  read_hdf5(file_id, "/YRECLDP_RCL_APB2", &yrecldp->rcl_apb2);
  read_hdf5(file_id, "/YRECLDP_RCL_APB3", &yrecldp->rcl_apb3);
  read_hdf5(file_id, "/YRECLDP_RCL_AS", &yrecldp->rcl_as);
  read_hdf5(file_id, "/YRECLDP_RCL_BS", &yrecldp->rcl_bs);
  read_hdf5(file_id, "/YRECLDP_RCL_CS", &yrecldp->rcl_cs);
  read_hdf5(file_id, "/YRECLDP_RCL_DS", &yrecldp->rcl_ds);
  read_hdf5(file_id, "/YRECLDP_RCL_X1S", &yrecldp->rcl_x1s);
  read_hdf5(file_id, "/YRECLDP_RCL_X2S", &yrecldp->rcl_x2s);
  read_hdf5(file_id, "/YRECLDP_RCL_X3S", &yrecldp->rcl_x3s);
  read_hdf5(file_id, "/YRECLDP_RCL_X4S", &yrecldp->rcl_x4s);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST1S", &yrecldp->rcl_const1s);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST2S", &yrecldp->rcl_const2s);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST3S", &yrecldp->rcl_const3s);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST4S", &yrecldp->rcl_const4s);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST5S", &yrecldp->rcl_const5s);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST6S", &yrecldp->rcl_const6s);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST7S", &yrecldp->rcl_const7s);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST8S", &yrecldp->rcl_const8s);
  read_hdf5(file_id, "/YRECLDP_RDENSWAT", &yrecldp->rdenswat);
  read_hdf5(file_id, "/YRECLDP_RDENSREF", &yrecldp->rdensref);
  read_hdf5(file_id, "/YRECLDP_RCL_AR", &yrecldp->rcl_ar);
  read_hdf5(file_id, "/YRECLDP_RCL_BR", &yrecldp->rcl_br);
  read_hdf5(file_id, "/YRECLDP_RCL_CR", &yrecldp->rcl_cr);
  read_hdf5(file_id, "/YRECLDP_RCL_DR", &yrecldp->rcl_dr);
  read_hdf5(file_id, "/YRECLDP_RCL_X1R", &yrecldp->rcl_x1r);
  read_hdf5(file_id, "/YRECLDP_RCL_X2R", &yrecldp->rcl_x2r);
  read_hdf5(file_id, "/YRECLDP_RCL_X4R", &yrecldp->rcl_x4r);
  read_hdf5(file_id, "/YRECLDP_RCL_KA273", &yrecldp->rcl_ka273);
  read_hdf5(file_id, "/YRECLDP_RCL_CDENOM1", &yrecldp->rcl_cdenom1);
  read_hdf5(file_id, "/YRECLDP_RCL_CDENOM2", &yrecldp->rcl_cdenom2);
  read_hdf5(file_id, "/YRECLDP_RCL_CDENOM3", &yrecldp->rcl_cdenom3);
  read_hdf5(file_id, "/YRECLDP_RCL_SCHMIDT", &yrecldp->rcl_schmidt);
  read_hdf5(file_id, "/YRECLDP_RCL_DYNVISC", &yrecldp->rcl_dynvisc);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST1R", &yrecldp->rcl_const1r);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST2R", &yrecldp->rcl_const2r);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST3R", &yrecldp->rcl_const3r);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST4R", &yrecldp->rcl_const4r);
  read_hdf5(file_id, "/YRECLDP_RCL_FAC1", &yrecldp->rcl_fac1);
  read_hdf5(file_id, "/YRECLDP_RCL_FAC2", &yrecldp->rcl_fac2);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST5R", &yrecldp->rcl_const5r);
  read_hdf5(file_id, "/YRECLDP_RCL_CONST6R", &yrecldp->rcl_const6r);
  read_hdf5(file_id, "/YRECLDP_RCL_FZRAB", &yrecldp->rcl_fzrab);
  read_hdf5(file_id, "/YRECLDP_RCL_FZRBB", &yrecldp->rcl_fzrbb);
  read_hdf5_int(file_id, "/YRECLDP_LCLDEXTRA", &yrecldp->lcldextra); // Bool
  read_hdf5_int(file_id, "/YRECLDP_LCLDBUDGET", &yrecldp->lcldbudget); // Bool
  read_hdf5_int(file_id, "/YRECLDP_NSSOPT", &yrecldp->nssopt);
  read_hdf5_int(file_id, "/YRECLDP_NCLDTOP", &yrecldp->ncldtop);
  read_hdf5_int(file_id, "/YRECLDP_NAECLBC", &yrecldp->naeclbc);
  read_hdf5_int(file_id, "/YRECLDP_NAECLDU", &yrecldp->naecldu);
  read_hdf5_int(file_id, "/YRECLDP_NAECLOM", &yrecldp->naeclom);
  read_hdf5_int(file_id, "/YRECLDP_NAECLSS", &yrecldp->naeclss);
  read_hdf5_int(file_id, "/YRECLDP_NAECLSU", &yrecldp->naeclsu);
  read_hdf5_int(file_id, "/YRECLDP_NCLDDIAG", &yrecldp->nclddiag);
  read_hdf5_int(file_id, "/YRECLDP_NAERCLD", &yrecldp->naercld);
  read_hdf5_int(file_id, "/YRECLDP_LAERLIQAUTOLSP", &yrecldp->laerliqautolsp); // Bool
  read_hdf5_int(file_id, "/YRECLDP_LAERLIQAUTOCP", &yrecldp->laerliqautocp); // Bool
  read_hdf5_int(file_id, "/YRECLDP_LAERLIQAUTOCPB", &yrecldp->laerliqautocpb); // Bool
  read_hdf5_int(file_id, "/YRECLDP_LAERLIQCOLL", &yrecldp->laerliqcoll); // Bool
  read_hdf5_int(file_id, "/YRECLDP_LAERICESED", &yrecldp->laericesed); // Bool
  read_hdf5_int(file_id, "/YRECLDP_LAERICEAUTO", &yrecldp->laericeauto); // Bool
  read_hdf5(file_id, "/YRECLDP_NSHAPEP", &yrecldp->nshapep);
  read_hdf5(file_id, "/YRECLDP_NSHAPEQ", &yrecldp->nshapeq);
  read_hdf5_int(file_id, "/YRECLDP_NBETA", &yrecldp->nbeta);

  status = H5Fclose(file_id);

#endif
}

/* Read reference result into memory */
void load_reference(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma,
		    dtype *plude, dtype *pcovptot, dtype *prainfrac_toprfz, dtype *pfsqlf, dtype *pfsqif,
		    dtype *pfcqlng, dtype *pfcqnng, dtype *pfsqrf, dtype *pfsqsf, dtype *pfcqrng, dtype *pfcqsng,
		    dtype *pfsqltur, dtype *pfsqitur, dtype *pfplsl, dtype *pfplsn, dtype *pfhpsl, dtype *pfhpsn,
		    dtype *tend_loc_a, dtype *tend_loc_q, dtype *tend_loc_t, dtype *tend_loc_cld)
{
  int nblocks = (ngptot / nproma) + min(ngptot % nproma, 1);

#ifdef HAVE_SERIALBOX
  serialboxSerializer_t* serializer = serialboxSerializerCreate(Read, "./data", "reference", "Binary");
  serialboxMetainfo_t* metainfo = serialboxSerializerGetGlobalMetainfo(serializer);
  serialboxSavepoint_t** savepoints = serialboxSerializerGetSavepointVector(serializer);
  serialboxSavepoint_t* savepoint = savepoints[0];

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
#endif

#ifdef HAVE_HDF5

  hid_t file_id, dataset_id;
  herr_t  status;
  file_id = H5Fopen(REFERENCE_FILE, H5F_ACC_RDONLY, H5P_DEFAULT);

  load_and_expand_2d(file_id, "PLUDE", nlon, nlev, nproma, ngptot, nblocks, plude);
  load_and_expand_2d(file_id, "PCOVPTOT", nlon, nlev, nproma, ngptot, nblocks, pcovptot);
  load_and_expand_1d(file_id, "PRAINFRAC_TOPRFZ", nlon, nproma, ngptot, nblocks, prainfrac_toprfz);
  load_and_expand_2d(file_id, "PFSQLF", nlon, nlev+1, nproma, ngptot, nblocks, pfsqlf);
  load_and_expand_2d(file_id, "PFSQIF", nlon, nlev+1, nproma, ngptot, nblocks, pfsqif);
  load_and_expand_2d(file_id, "PFCQLNG", nlon, nlev+1, nproma, ngptot, nblocks, pfcqlng);
  load_and_expand_2d(file_id, "PFCQNNG", nlon, nlev+1, nproma, ngptot, nblocks, pfcqnng);
  load_and_expand_2d(file_id, "PFSQRF", nlon, nlev+1, nproma, ngptot, nblocks, pfsqrf);
  load_and_expand_2d(file_id, "PFSQSF", nlon, nlev+1, nproma, ngptot, nblocks, pfsqsf);
  load_and_expand_2d(file_id, "PFCQRNG", nlon, nlev+1, nproma, ngptot, nblocks, pfcqrng);
  load_and_expand_2d(file_id, "PFCQSNG", nlon, nlev+1, nproma, ngptot, nblocks, pfcqsng);
  load_and_expand_2d(file_id, "PFSQLTUR", nlon, nlev+1, nproma, ngptot, nblocks, pfsqltur);
  load_and_expand_2d(file_id, "PFSQITUR", nlon, nlev+1, nproma, ngptot, nblocks, pfsqitur);
  load_and_expand_2d(file_id, "PFPLSL", nlon, nlev+1, nproma, ngptot, nblocks, pfplsl);
  load_and_expand_2d(file_id, "PFPLSN", nlon, nlev+1, nproma, ngptot, nblocks, pfplsn);
  load_and_expand_2d(file_id, "PFHPSL", nlon, nlev+1, nproma, ngptot, nblocks, pfhpsl);
  load_and_expand_2d(file_id, "PFHPSN", nlon, nlev+1, nproma, ngptot, nblocks, pfhpsn);
  load_and_expand_2d(file_id, "TENDENCY_LOC_T", nlon, nlev, nproma, ngptot, nblocks, tend_loc_t);
  load_and_expand_2d(file_id, "TENDENCY_LOC_Q", nlon, nlev, nproma, ngptot, nblocks, tend_loc_q);
  load_and_expand_2d(file_id, "TENDENCY_LOC_A", nlon, nlev, nproma, ngptot, nblocks, tend_loc_a);
  load_and_expand_3d(file_id, "TENDENCY_LOC_CLD", nlon, nlev, nclv, nproma, ngptot, nblocks, tend_loc_cld);

  status = H5Fclose(file_id);
#endif
}
