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

void read_hdf5(hid_t file_id, const char *name, double *field) {
  hid_t dataset_id;
  herr_t  status;
  dataset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, field);
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
  file_id = H5Fopen(INPUT_FILE, H5F_ACC_RDWR, H5P_DEFAULT);
  
  read_hdf5_int(file_id, "/KLEV", klev);
  read_hdf5_int(file_id, "/KLON", klon);
  
  status = H5Fclose(file_id);
#endif
}

void expand_1d(double *buffer, double *field_in, int nlon, int nproma, int ngptot, int nblocks)
{
  int b, l, i, buf_start_idx, buf_idx;
  double (*field)[nproma] = (double (*)[nproma]) field_in;

#pragma omp parallel for default(shared) private(b, l, i, buf_start_idx, buf_idx)
  for (b = 0; b < nblocks; b++) {
    buf_start_idx = ((b)*nproma) % nlon;
    for (i = 0; i < nproma; i++) {
      buf_idx = (buf_start_idx + i) % nlon;
      field[b][i] = buffer[buf_idx];
    }
  }

  // Zero out the remainder of last block
  /*
  bsize = min(nproma, ngptot - (nblocks-1)*nproma);  // Size of the field block
  printf("zeroing last block : %d \n",bsize);
  for (i=bsize; i<nproma; i++) {
    for (l=0; l<nlev; l++) {
      field[nblocks-1][l][i] = 0.0;
    }
  }
  */

}


void expand_1d_int(int *buffer, int *field_in, int nlon, int nproma, int ngptot, int nblocks)
{
  int b, l, i, buf_start_idx, buf_idx;
  int (*field)[nproma] = (int (*)[nproma]) field_in;

  #pragma omp parallel for default(shared) private(b, l, i, buf_start_idx, buf_idx)
  for (b = 0; b < nblocks; b++) {
    buf_start_idx = ((b)*nproma) % nlon;
    for (i = 0; i < nproma; i++) {
      buf_idx = (buf_start_idx + i) % nlon;
      field[b][i] = buffer[buf_idx];
    }
  }

  // Zero out the remainder of last block
  /*
  bsize = min(nproma, ngptot - (nblocks-1)*nproma);  // Size of the field block
  printf("zeroing last block : %d \n",bsize);
  for (i=bsize; i<nproma; i++) {
    for (l=0; l<nlev; l++) {
      field[nblocks-1][l][i] = 0;
    }
  }
  */
}

void expand_2d(double *buffer_in, double *field_in, int nlon, int nlev, int nproma, int ngptot, int nblocks)
{
  int b, l, i, buf_start_idx, buf_idx;
  double (*buffer)[nlon] = (double (*)[nlon]) buffer_in;
  double (*field)[nlev][nproma] = (double (*)[nlev][nproma]) field_in;
  
  #pragma omp parallel for default(shared) private(b, buf_start_idx, buf_idx, l, i)
  for (b = 0; b < nblocks; b++) {
    buf_start_idx = ((b)*nproma) % nlon;
    for (i = 0; i < nproma; i++) {
    for (l = 0; l < nlev; l++) {
       buf_idx = (buf_start_idx + i) % nlon;
       field[b][l][i] = buffer[l][buf_idx];       
    }
    }
  }

  // Zero out the remainder of last block
  /*
  bsize = min(nproma, ngptot - (nblocks-1)*nproma);  // Size of the field block
  printf("zeroing last block : %d \n",bsize);
  for (i=bsize; i<nproma; i++) {
    for (l=0; l<nlev; l++) {
      field[nblocks-1][l][i] = 0.0;
    }
  }
  */
  
}


void expand_3d(double *buffer_in, double *field_in, int nlon, int nlev, int nclv, int nproma, int ngptot, int nblocks)
{
  int b, l, c, i, buf_start_idx, buf_idx;
  double (*buffer)[nlev][nlon] = (double (*)[nlev][nlon]) buffer_in;
  double (*field)[nclv][nlev][nproma] = (double (*)[nclv][nlev][nproma]) field_in;
  
#pragma omp parallel for default(shared) private(b, buf_start_idx, buf_idx, l, i)
  for (b = 0; b < nblocks; b++) {
    buf_start_idx = ((b)*nproma) % nlon;   
    for (i = 0; i < nproma; i++) {
    for (c = 0; c < nclv; c++) {
    for (l = 0; l < nlev; l++) {
      buf_idx = (buf_start_idx + i) % nlon;
      field[b][c][l][i] = buffer[c][l][buf_idx];
    }
    }
    }
  }

  // Zero out the remainder of last block
  /*
  bsize = min(nproma, ngptot - (nblocks-1)*nproma);  // Size of the field block
  printf("zeroing last block : %d \n",bsize);
  for (i=bsize; i<nproma; i++) {
    for (l=0; l<nlev; l++) {
      field[nblocks-1][l][i] = 0.0;
    }
  }
  */

}

#ifdef HAVE_SERIALBOX
void load_and_expand_1d(serialboxSerializer_t *serializer, serialboxSavepoint_t* savepoint,
    const char *name, int nlon, int nproma, int ngptot, int nblocks, double *field)
{
  double buffer[nlon];
  int strides[1] = {1};

  serialboxSerializerRead(serializer, name, savepoint, buffer, strides, 1);
  expand_1d((double *)buffer, field, nlon, nproma, ngptot, nblocks);
}

void load_and_expand_1d_int(serialboxSerializer_t *serializer, serialboxSavepoint_t* savepoint,
    const char *name, int nlon, int nproma, int ngptot, int nblocks, int *field)
{
  int buffer[nlon];
  int strides[1] = {1};

  serialboxSerializerRead(serializer, name, savepoint, buffer, strides, 1);
  expand_1d_int((int *)buffer, field, nlon, nproma, ngptot, nblocks);
}

void load_and_expand_2d(serialboxSerializer_t *serializer, serialboxSavepoint_t* savepoint,
    const char *name, int nlon, int nlev, int nproma, int ngptot, int nblocks, double *field)
{
  double buffer[nlev][nlon];
  int strides[2] = {1, nlon};

  serialboxSerializerRead(serializer, name, savepoint, buffer, strides, 2);
  expand_2d((double *)buffer, field, nlon, nlev, nproma, ngptot, nblocks);
}

void load_and_expand_3d(serialboxSerializer_t *serializer, serialboxSavepoint_t* savepoint,
    const char *name, int nlon, int nlev, int nclv, int nproma, int ngptot, int nblocks, double *field)
{
  double buffer[nclv][nlev][nlon];
  int strides[3] = {1, nlon, nlev*nlon};

  serialboxSerializerRead(serializer, name, savepoint, buffer, strides, 3);
  expand_3d((double *)buffer, field, nlon, nlev, nclv, nproma, ngptot, nblocks);
}
#endif

#if HAVE_HDF5
void load_and_expand_1d(hid_t file_id, const char *name, int nlon, int nproma, int ngptot, int nblocks, double *field)
{
  double buffer[nlon];
  int strides[1] = {1};
  hid_t dataset_id;
  dataset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
  herr_t  status;
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  status = H5Dclose(dataset_id);
  expand_1d((double *)buffer, field, nlon, nproma, ngptot, nblocks);
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

void load_and_expand_2d(hid_t file_id, const char *name, int nlon, int nlev, int nproma, int ngptot, int nblocks, double *field)
{
  double buffer[nlev][nlon];
  int strides[2] = {1, nlon};
  hid_t dataset_id;
  dataset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
  herr_t  status;
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  status = H5Dclose(dataset_id);
  expand_2d((double *)buffer, field, nlon, nlev, nproma, ngptot, nblocks);
}

void load_and_expand_3d(hid_t file_id, const char *name, int nlon, int nlev, int nclv, int nproma, int ngptot, int nblocks, double *field)
{
  double buffer[nclv][nlev][nlon];
  int strides[3] = {1, nlon, nlev*nlon};
  hid_t dataset_id;
  dataset_id = H5Dopen2(file_id, name, H5P_DEFAULT);
  herr_t  status;
  status = H5Dread(dataset_id, H5T_NATIVE_DOUBLE, H5S_ALL, H5S_ALL, H5P_DEFAULT, buffer);
  status = H5Dclose(dataset_id);
  expand_3d((double *)buffer, field, nlon, nlev, nclv, nproma, ngptot, nblocks);
}
#endif

/* Read input state into memory */
void load_state(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma, 
    double* ptsphy, double* plcrit_aer, double* picrit_aer,
    double* pre_ice, double* pccn, double* pnice, double* pt, double* pq, 
    double* tend_cml_t, double* tend_cml_q, double* tend_cml_a, double* tend_cml_cld,
    double* tend_tmp_t, double* tend_tmp_q, double* tend_tmp_a, double* tend_tmp_cld,
    double* pvfa, double* pvfl, double* pvfi, double* pdyna, double* pdynl, double* pdyni, 
    double* phrsw, double* phrlw, double* pvervel, double* pap, double* paph, double* plsm,
    int* ktype, double* plu, double* plude, double* psnde, double* pmfu,
    double* pmfd, double* pa, double* pclv, double* psupsat)
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

  *ptsphy = serialboxMetainfoGetFloat64(metainfo, "PTSPHY");

  /* Populate global parameter values from meta-data */
  rg = serialboxMetainfoGetFloat64(metainfo, "RG");
  rd = serialboxMetainfoGetFloat64(metainfo, "RD");
  rcpd = serialboxMetainfoGetFloat64(metainfo, "RCPD");
  retv = serialboxMetainfoGetFloat64(metainfo, "RETV");
  rlvtt = serialboxMetainfoGetFloat64(metainfo, "RLVTT");
  rlstt = serialboxMetainfoGetFloat64(metainfo, "RLSTT");
  rlmlt = serialboxMetainfoGetFloat64(metainfo, "RLMLT");
  rtt = serialboxMetainfoGetFloat64(metainfo, "RTT");
  rv = serialboxMetainfoGetFloat64(metainfo, "RV");
  r2es = serialboxMetainfoGetFloat64(metainfo, "R2ES");
  r3les = serialboxMetainfoGetFloat64(metainfo, "R3LES");
  r3ies = serialboxMetainfoGetFloat64(metainfo, "R3IES");
  r4les = serialboxMetainfoGetFloat64(metainfo, "R4LES");
  r4ies = serialboxMetainfoGetFloat64(metainfo, "R4IES");
  r5les = serialboxMetainfoGetFloat64(metainfo, "R5LES");
  r5ies = serialboxMetainfoGetFloat64(metainfo, "R5IES");
  r5alvcp = serialboxMetainfoGetFloat64(metainfo, "R5ALVCP");
  r5alscp = serialboxMetainfoGetFloat64(metainfo, "R5ALSCP");
  ralvdcp = serialboxMetainfoGetFloat64(metainfo, "RALVDCP");
  ralsdcp = serialboxMetainfoGetFloat64(metainfo, "RALSDCP");
  ralfdcp = serialboxMetainfoGetFloat64(metainfo, "RALFDCP");
  rtwat = serialboxMetainfoGetFloat64(metainfo, "RTWAT");
  rtice = serialboxMetainfoGetFloat64(metainfo, "RTICE");
  rticecu = serialboxMetainfoGetFloat64(metainfo, "RTICECU");
  rtwat_rtice_r = serialboxMetainfoGetFloat64(metainfo, "RTWAT_RTICE_R");
  rtwat_rticecu_r = serialboxMetainfoGetFloat64(metainfo, "RTWAT_RTICECU_R");
  rkoop1 = serialboxMetainfoGetFloat64(metainfo, "RKOOP1");
  rkoop2 = serialboxMetainfoGetFloat64(metainfo, "RKOOP2");

  yrecldp->ramid = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RAMID");
  yrecldp->rcldiff = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLDIFF");
  yrecldp->rcldiff_convi = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLDIFF_CONVI");
  yrecldp->rclcrit = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLCRIT");
  yrecldp->rclcrit_sea = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLCRIT_SEA");
  yrecldp->rclcrit_land = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLCRIT_LAND");
  yrecldp->rkconv = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RKCONV");
  yrecldp->rprc1 = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RPRC1");
  yrecldp->rprc2 = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RPRC2");
  yrecldp->rcldmax = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLDMAX");
  yrecldp->rpecons = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RPECONS");
  yrecldp->rvrfactor = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RVRFACTOR");
  yrecldp->rprecrhmax = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RPRECRHMAX");
  yrecldp->rtaumel = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RTAUMEL");
  yrecldp->ramin = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RAMIN");
  yrecldp->rlmin = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RLMIN");
  yrecldp->rkooptau = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RKOOPTAU");

  yrecldp->rcldtopp = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLDTOPP");
  yrecldp->rlcritsnow = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RLCRITSNOW");
  yrecldp->rsnowlin1 = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RSNOWLIN1");
  yrecldp->rsnowlin2 = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RSNOWLIN2");
  yrecldp->ricehi1 = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RICEHI1");
  yrecldp->ricehi2 = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RICEHI2");
  yrecldp->riceinit = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RICEINIT");
  yrecldp->rvice = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RVICE");
  yrecldp->rvrain = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RVRAIN");
  yrecldp->rvsnow = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RVSNOW");
  yrecldp->rthomo = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RTHOMO");
  yrecldp->rcovpmin = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCOVPMIN");
  yrecldp->rccn = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCCN");
  yrecldp->rnice = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RNICE");
  yrecldp->rccnom = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCCNOM");
  yrecldp->rccnss = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCCNSS");
  yrecldp->rccnsu = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCCNSU");
  yrecldp->rcldtopcf = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCLDTOPCF");
  yrecldp->rdepliqrefrate = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RDEPLIQREFRATE");
  yrecldp->rdepliqrefdepth = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RDEPLIQREFDEPTH");
  yrecldp->rcl_kkaac = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KKAac");
  yrecldp->rcl_kkbac = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KKBac");
  yrecldp->rcl_kkaau = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KKAau");
  yrecldp->rcl_kkbauq = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KKBauq");
  yrecldp->rcl_kkbaun = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KKBaun");
  yrecldp->rcl_kk_cloud_num_sea = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KK_cloud_num_sea");
  yrecldp->rcl_kk_cloud_num_land = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KK_cloud_num_land");
  yrecldp->rcl_ai = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_AI");
  yrecldp->rcl_bi = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_BI");
  yrecldp->rcl_ci = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CI");
  yrecldp->rcl_di = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_DI");
  yrecldp->rcl_x1i = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X1I");
  yrecldp->rcl_x2i = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X2I");
  yrecldp->rcl_x3i = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X3I");
  yrecldp->rcl_x4i = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X4I");
  yrecldp->rcl_const1i = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST1I");
  yrecldp->rcl_const2i = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST2I");
  yrecldp->rcl_const3i = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST3I");
  yrecldp->rcl_const4i = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST4I");
  yrecldp->rcl_const5i = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST5I");
  yrecldp->rcl_const6i = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST6I");
  yrecldp->rcl_apb1 = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_APB1");
  yrecldp->rcl_apb2 = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_APB2");
  yrecldp->rcl_apb3 = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_APB3");
  yrecldp->rcl_as = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_AS");
  yrecldp->rcl_bs = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_BS");
  yrecldp->rcl_cs = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CS");
  yrecldp->rcl_ds = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_DS");
  yrecldp->rcl_x1s = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X1S");
  yrecldp->rcl_x2s = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X2S");
  yrecldp->rcl_x3s = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X3S");
  yrecldp->rcl_x4s = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X4S");
  yrecldp->rcl_const1s = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST1S");
  yrecldp->rcl_const2s = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST2S");
  yrecldp->rcl_const3s = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST3S");
  yrecldp->rcl_const4s = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST4S");
  yrecldp->rcl_const5s = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST5S");
  yrecldp->rcl_const6s = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST6S");
  yrecldp->rcl_const7s = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST7S");
  yrecldp->rcl_const8s = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST8S");
  yrecldp->rdenswat = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RDENSWAT");
  yrecldp->rdensref = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RDENSREF");
  yrecldp->rcl_ar = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_AR");
  yrecldp->rcl_br = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_BR");
  yrecldp->rcl_cr = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CR");
  yrecldp->rcl_dr = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_DR");
  yrecldp->rcl_x1r = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X1R");
  yrecldp->rcl_x2r = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X2R");
  yrecldp->rcl_x4r = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_X4R");
  yrecldp->rcl_ka273 = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_KA273");
  yrecldp->rcl_cdenom1 = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CDENOM1");
  yrecldp->rcl_cdenom2 = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CDENOM2");
  yrecldp->rcl_cdenom3 = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CDENOM3");
  yrecldp->rcl_schmidt = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_SCHMIDT");
  yrecldp->rcl_dynvisc = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_DYNVISC");
  yrecldp->rcl_const1r = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST1R");
  yrecldp->rcl_const2r = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST2R");
  yrecldp->rcl_const3r = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST3R");
  yrecldp->rcl_const4r = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST4R");
  yrecldp->rcl_fac1 = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_FAC1");
  yrecldp->rcl_fac2 = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_FAC2");
  yrecldp->rcl_const5r = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST5R");
  yrecldp->rcl_const6r = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_CONST6R");
  yrecldp->rcl_fzrab = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_FZRAB");
  yrecldp->rcl_fzrbb = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_RCL_FZRBB");
  yrecldp->lcldextra = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LCLDEXTRA");
  yrecldp->lcldbudget = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LCLDBUDGET");
  yrecldp->nssopt = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NSSOPT");
  yrecldp->ncldtop = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NCLDTOP");
  yrecldp->naeclbc = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAECLBC");
  yrecldp->naecldu = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAECLDU");
  yrecldp->naeclom = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAECLOM");
  yrecldp->naeclss = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAECLSS");
  yrecldp->naeclsu = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAECLSU");
  yrecldp->nclddiag = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NCLDDIAG");
  yrecldp->naercld = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NAERCLD");
  yrecldp->laerliqautolsp = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERLIQAUTOLSP");
  yrecldp->laerliqautocp = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERLIQAUTOCP");
  yrecldp->laerliqautocpb = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERLIQAUTOCPB");
  yrecldp->laerliqcoll = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERLIQCOLL");
  yrecldp->laericesed = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERICESED");
  yrecldp->laericeauto = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_LAERICEAUTO");
  yrecldp->nshapep = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NSHAPEP");
  yrecldp->nshapeq = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NSHAPEQ");
  yrecldp->nbeta = serialboxMetainfoGetFloat64(metainfo, "YRECLDP_NBETA");

  serialboxSerializerDestroySavepointVector(savepoints, 1);
  serialboxSerializerDestroy(serializer);
#endif

#ifdef HAVE_HDF5

  hid_t file_id, dataset_id;
  herr_t  status;
  file_id = H5Fopen(INPUT_FILE, H5F_ACC_RDWR, H5P_DEFAULT);

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
  
  read_hdf5(file_id, "/RG", &rg);
  read_hdf5(file_id, "/RD", &rd);
  read_hdf5(file_id, "/RCPD", &rcpd); 
  read_hdf5(file_id, "/RETV", &retv);   
  read_hdf5(file_id, "/RLVTT", &rlvtt);  
  read_hdf5(file_id, "/RLSTT", &rlstt);  
  read_hdf5(file_id, "/RLMLT", &rlmlt);  
  read_hdf5(file_id, "/RTT", &rtt);  
  read_hdf5(file_id, "/RV", &rv);  
  read_hdf5(file_id, "/R2ES", &r2es);  
  read_hdf5(file_id, "/R3LES", &r3les);  
  read_hdf5(file_id, "/R3IES", &r3ies);  
  read_hdf5(file_id, "/R4LES", &r4les);  
  read_hdf5(file_id, "/R4IES", &r4ies);  
  read_hdf5(file_id, "/R5LES", &r5les);  
  read_hdf5(file_id, "/R5IES", &r5ies);  
  read_hdf5(file_id, "/R5ALVCP", &r5alvcp);  
  read_hdf5(file_id, "/R5ALSCP", &r5alscp);  
  read_hdf5(file_id, "/RALVDCP", &ralvdcp);  
  read_hdf5(file_id, "/RALSDCP", &ralsdcp);  
  read_hdf5(file_id, "/RALFDCP", &ralfdcp);  
  read_hdf5(file_id, "/RTWAT", &rtwat);  
  read_hdf5(file_id, "/RTICE", &rtice);  
  read_hdf5(file_id, "/RTICECU", &rticecu);  
  read_hdf5(file_id, "/RTWAT_RTICE_R", &rtwat_rtice_r);  
  read_hdf5(file_id, "/RTWAT_RTICECU_R", &rtwat_rticecu_r);  
  read_hdf5(file_id, "/RKOOP1", &rkoop1);  
  read_hdf5(file_id, "/RKOOP2", &rkoop2);  
  
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
		    double *plude, double *pcovptot, double *prainfrac_toprfz, double *pfsqlf, double *pfsqif,
		    double *pfcqlng, double *pfcqnng, double *pfsqrf, double *pfsqsf, double *pfcqrng, double *pfcqsng,
		    double *pfsqltur, double *pfsqitur, double *pfplsl, double *pfplsn, double *pfhpsl, double *pfhpsn,
		    double *tend_loc_a, double *tend_loc_q, double *tend_loc_t, double *tend_loc_cld)
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
  file_id = H5Fopen(REFERENCE_FILE, H5F_ACC_RDWR, H5P_DEFAULT);

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
