#include "load_state.h"

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

void expand_1d(double *buffer, double *field_in, int nlon, int nproma, int ngptot, int nblocks)
{
  int gidx, bsize, bidx, bend, fsize, b, l, i;
  double (*field)[nproma] = (double (*)[nproma]) field_in;

  #pragma omp parallel for default(shared) private(gidx, bsize, bidx, bend, fsize, b, l, i)
  for (b = 0; b < nblocks; b++) {
    gidx = b*nproma;  // Global starting index of the block in the general domain
    bsize = min(nproma, ngptot - gidx);  // Size of the field block
    bidx = gidx % nlon;  // Rolling index into the input buffer
    bend = min(bidx+bsize, nlon);  // Idealised final index in the input buffer

    if (bend-bidx < bsize) {
      // The input buffer does not hold enough data to fill field block;
      // we need to fill the rest of the block with data from front of buffer.
      fsize = nlon - bidx;
      for (i = 0; i < fsize; i++) { field[b][i] = buffer[bidx + i]; }
      for (i = 0; i < bsize-fsize; i++) { field[b][fsize+i] = buffer[i]; }
    } else {
      // Simply copy a block of data from the rolling buffer index
      for (i = 0; i < bsize; i++) { field[b][i] = buffer[bidx+i]; }
    }
    // Zero out the remainder of last block
    for (i=bsize; i<nproma; i++) { field[b][i] = 0.0; }
  }
}

void expand_1d_int(int *buffer, int *field_in, int nlon, int nproma, int ngptot, int nblocks)
{
  int gidx, bsize, bidx, bend, fsize, b, l, i;
  int (*field)[nproma] = (int (*)[nproma]) field_in;

  #pragma omp parallel for default(shared) private(gidx, bsize, bidx, bend, fsize, b, l, i)
  for (b = 0; b < nblocks; b++) {
    gidx = b*nproma;  // Global starting index of the block in the general domain
    bsize = min(nproma, ngptot - gidx);  // Size of the field block
    bidx = gidx % nlon;  // Rolling index into the input buffer
    bend = min(bidx+bsize, nlon);  // Idealised final index in the input buffer

    if (bend-bidx < bsize) {
      // The input buffer does not hold enough data to fill field block;
      // we need to fill the rest of the block with data from front of buffer.
      fsize = nlon - bidx;
      for (i = 0; i < fsize; i++) { field[b][i] = buffer[bidx + i]; }
      for (i = 0; i < bsize-fsize; i++) { field[b][fsize+i] = buffer[i]; }
    } else {
      // Simply copy a block of data from the rolling buffer index
      for (i = 0; i < bsize; i++) { field[b][i] = buffer[bidx+i]; }
    }
    // Zero out the remainder of last block
    for (i=bsize; i<nproma; i++) { field[b][i] = 0.0; }
  }
}

void expand_2d(double *buffer_in, double *field_in, int nlon, int nlev, int nproma, int ngptot, int nblocks)
{
  int gidx, bsize, bidx, bend, fsize, b, l, i;
  double (*buffer)[nlon] = (double (*)[nlon]) buffer_in;
  double (*field)[nlev][nproma] = (double (*)[nlev][nproma]) field_in;

  #pragma omp parallel for default(shared) private(gidx, bsize, bidx, bend, fsize, b, l, i)
  for (b = 0; b < nblocks; b++) {
    gidx = b*nproma;  // Global starting index of the block in the general domain
    bsize = min(nproma, ngptot - gidx);  // Size of the field block
    bidx = gidx % nlon;  // Rolling index into the input buffer
    bend = min(bidx+bsize, nlon);  // Idealised final index in the input buffer

    if (bend-bidx < bsize) {
      // The input buffer does not hold enough data to fill field block;
      // we need to fill the rest of the block with data from front of buffer.
      fsize = nlon - bidx;
      for (l = 0; l < nlev; l++) {
      	for (i = 0; i < fsize; i++) { field[b][l][i] = buffer[l][bidx + i]; }
	for (i = 0; i < bsize-fsize; i++) { field[b][l][fsize+i] = buffer[l][i]; }
      }
    } else {
      // Simply copy a block of data from the rolling buffer index
      for (l = 0; l < nlev; l++) {
	for (i = 0; i < bsize; i++) { field[b][l][i] = buffer[l][bidx+i]; }
      }
    }
    // Zero out the remainder of last block
    for (l=0; l<nlev; l++ ) { 
      for (i=bsize; i<nproma; i++) { field[b][l][i] = 0.0; }
    }
  }
}

void expand_3d(double *buffer_in, double *field_in, int nlon, int nlev, int nclv, int nproma, int ngptot, int nblocks)
{
  int gidx, bsize, bidx, bend, fsize, b, l, c, i;
  double (*buffer)[nlev][nlon] = (double (*)[nlev][nlon]) buffer_in;
  double (*field)[nclv][nlev][nproma] = (double (*)[nclv][nlev][nproma]) field_in;

#pragma omp parallel for default(shared) private(gidx, bsize, bidx, bend, fsize, b, l, c, i)
  for (b = 0; b < nblocks; b++) {
    gidx = b*nproma;  // Global starting index of the block in the general domain
    bsize = min(nproma, ngptot - gidx);  // Size of the field block
    bidx = gidx % nlon;  // Rolling index into the input buffer
    bend = min(bidx+bsize, nlon);  // Idealised final index in the input buffer

    if (bend-bidx < bsize) {
      // The input buffer does not hold enough data to fill field block;
      // we need to fill the rest of the block with data from front of buffer.
      fsize = nlon - bidx;
      for (c = 0; c < nclv; c++) {
	for (l = 0; l < nlev; l++) {
	  for (i = 0; i < fsize; i++) { field[b][c][l][i] = buffer[c][l][bidx + i]; }
	  for (i = 0; i < bsize-fsize; i++) { field[b][c][l][fsize+i] = buffer[c][l][i]; }
	}
      }
    } else {
      // Simply copy a block of data from the rolling buffer index
      for (c = 0; c < nclv; c++) {
	for (l = 0; l < nlev; l++) {
	  for (i = 0; i < bsize; i++) { field[b][c][l][i] = buffer[c][l][bidx+i]; }
	}
      }
    }
    // Zero out the remainder of last block
    for (c = 0; c < nclv; c++) {
      for (l=0; l < nlev; l++ ) { 
	for (i=bsize; i<nproma; i++) { field[b][c][l][i] = 0.0; }
      }
    }
  }
}


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


/* Read input state into memory */
void load_state(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma, 
    double* ptsphy, double* plcrit_aer, double* picrit_aer,
    double* pre_ice, double* pccn, double* pnice, double* pt, double* pq, 
    double* tend_cml_t, double* tend_cml_q, double* tend_cml_a, double* tend_cml_cld,
    double* tend_tmp_t, double* tend_tmp_q, double* tend_tmp_a, double* tend_tmp_cld,
    double* pvfa, double* pvfl, double* pvfi, double* pdyna, double* pdynl, double* pdyni, 
    double* phrsw, double* phrlw, double* pvervel, double* pap, double* paph, double* plsm,
    int* ldcum, int* ktype, double* plu, double* plude, double* psnde, double* pmfu,
    double* pmfd, double* pa, double* pclv, double* psupsat)
{
  serialboxSerializer_t* serializer = serialboxSerializerCreate(Read, "./data", "input", "Binary");
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
  load_and_expand_1d_int(serializer, savepoint, "LDCUM", nlon, nproma, ngptot, nblocks, ldcum);
  load_and_expand_1d_int(serializer, savepoint, "KTYPE", nlon, nproma, ngptot, nblocks, ktype);
  load_and_expand_2d(serializer, savepoint, "PLU", nlon, nlev, nproma, ngptot, nblocks, plu);
  load_and_expand_2d(serializer, savepoint, "PLUDE", nlon, nlev, nproma, ngptot, nblocks, plude);
  load_and_expand_2d(serializer, savepoint, "PSNDE", nlon, nlev, nproma, ngptot, nblocks, psnde);
  load_and_expand_2d(serializer, savepoint, "PMFU", nlon, nlev, nproma, ngptot, nblocks, pmfu);
  load_and_expand_2d(serializer, savepoint, "PMFD", nlon, nlev, nproma, ngptot, nblocks, pmfd);
  load_and_expand_2d(serializer, savepoint, "PA", nlon, nlev, nproma, ngptot, nblocks, pa);
  load_and_expand_3d(serializer, savepoint, "PCLV", nlon, nlev, nclv, nproma, ngptot, nblocks, pclv);
  load_and_expand_2d(serializer, savepoint, "PSUPSAT", nlon, nlev, nproma, ngptot, nblocks, psupsat);

  serialboxSerializerDestroySavepointVector(savepoints, 1);
  serialboxSerializerDestroy(serializer);
}

/*
 * Retrieve parametersand timestep size from the serializer
 */
void initialise_parameters(double ptsphy)
{
}
