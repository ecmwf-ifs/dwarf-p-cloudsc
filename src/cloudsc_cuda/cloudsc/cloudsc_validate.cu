/*
 * (C) Copyright 1988- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include "cloudsc_validate.h"

#include <float.h>
#include <math.h>

#define min(a, b) (((a) < (b)) ? (a) : (b))
#define max(a, b) (((a) > (b)) ? (a) : (b))


void print_error(const char *name, double zminval, double zmaxval, double zmaxerr,
		 double zerrsum, double zsum, double zavgpgp, int ndim)
{
  double zrelerr, zeps = DBL_EPSILON;
  int iopt = 0;
  if (zerrsum < zeps) {
    zrelerr = 0.0;
    iopt = 1;
  } else if (zsum < zeps) {
    zrelerr = zerrsum / (1.0 + zsum);
    iopt = 2;
  } else {
    zrelerr = zerrsum / zsum;
    iopt = 3;
  }

  //-- If you get 4 exclamation marks next to your error output,
  //   then it is likely that some uninitialized variables exists or
  //   some other screw-up -- watch out this !!!!
  //char *clwarn;
  const char* clwarn = (zrelerr > 10.0 * zeps) ? " !!!!" : "     ";
  zrelerr = 100.0 * zrelerr;

  printf(" zsum: %20.13le   zerrsum: %20.13le", zsum, zerrsum);
  printf(" %+20s %dD%d %20.13le %20.13le %20.13le %20.13le %20.13le %s\n",
         name, ndim, iopt, zminval, zmaxval, zmaxerr, zavgpgp, zrelerr, clwarn);
}


void validate_1d(const char *name, float * v_ref, float * v_field, int nlon, int ngptot, int nblocks)
{
  /* Computes and prints errors in the "L2 norm sense" */
  int b, bsize, jk;
  double zminval, zmaxval, zdiff, zmaxerr, zerrsum, zsum, zrelerr, zavgpgp;
  double field_val, ref_val;
  float (*field)[nlon] = (float (*)[nlon]) v_field;
  float (*reference)[nlon] = (float (*)[nlon]) v_ref;

  zminval = +DBL_MAX;
  zmaxval = -DBL_MAX;
  zmaxerr = 0.0;
  zerrsum = 0.0;
  zsum = 0.0;

  // # pragma omp parallel for default(shared) private(b, bsize, jk)		\
  //   reduction(min:zminval) reduction(max:zmaxval,zmaxerr) reduction(+:zerrsum,zsum)
  for (b = 0; b < nblocks; b++) {
    bsize = min(nlon, ngptot - b*nlon);  // field block size
    for (jk = 0; jk < bsize; jk++) {
      field_val = (double) field[b][jk];
      ref_val = (double) reference[b][jk];

      zminval = fmin(zminval, field_val);
      zmaxval = fmax(zmaxval, field_val);

      // Difference against reference result in one-norm sense
      zdiff = fabs(field_val - ref_val);
      zmaxerr = fmax(zmaxerr, zdiff);
      zerrsum = zerrsum + zdiff;
      zsum = zsum + fabs(ref_val);
    }
  }
  zavgpgp = zerrsum / (double) ngptot;
  print_error(name, zminval, zmaxval, zmaxerr, zerrsum, zsum, zavgpgp, 2);
}


void validate_2d(const char *name, float *v_ref, float *v_field, int nlon, int nlev, int ngptot, int nblocks)
{
  /* Computes and prints errors in the "L2 norm sense" */
  int b, bsize, jl, jk;
  double zminval, zmaxval, zdiff, zmaxerr, zerrsum, zsum, zrelerr, zavgpgp, zdiffmax;
  double field_val, ref_val;
  float (*field)[nlev][nlon] = (float (*)[nlev][nlon]) v_field;
  float (*reference)[nlev][nlon] = (float (*)[nlev][nlon]) v_ref;

  zminval = +DBL_MAX;
  zmaxval = -DBL_MAX;
  zmaxerr = 0.0;
  zerrsum = 0.0;
  zsum = 0.0;
  zdiffmax = 0.0;

  // # pragma omp parallel for default(shared) private(b, bsize, jl, jk) \
  //   reduction(min:zminval) reduction(max:zmaxval,zmaxerr) reduction(+:zerrsum,zsum)
  for (b = 0; b < nblocks; b++) {
    bsize = min(nlon, ngptot - b*nlon);  // field block size
    for (jl = 0; jl < nlev; jl++) {
      for (jk = 0; jk < bsize; jk++) {
        field_val = (double) field[b][jl][jk];
        ref_val = (double) reference[b][jl][jk];

	zminval = fmin(zminval, field_val);
	zmaxval = fmax(zmaxval, field_val);
	// Difference against reference result in one-norm sense
	zdiff = fabs(field_val - ref_val);
	zmaxerr = fmax(zmaxerr, zdiff);
	zdiffmax = fmax(zdiffmax, zdiff);
	zerrsum = zerrsum + zdiff;
	zsum = zsum + fabs(ref_val);
      }
    }
  }
  zavgpgp = zerrsum / (double) ngptot;
  print_error(name, zminval, zmaxval, zmaxerr, zerrsum, zsum, zavgpgp, 2);
}


void validate_3d(const char *name, float *v_ref, float *v_field, int nlon,
    int nlev, int nclv, int ngptot, int nblocks)
{
  /* Computes and prints errors in the "L2 norm sense" */
  int b, bsize, jl, jk, jm;
  double zminval, zmaxval, zdiff, zmaxerr, zerrsum, zsum, zrelerr, zavgpgp, zdiffmax;
  double field_val, ref_val;
  float (*field)[nclv][nlev][nlon] = (float (*)[nclv][nlev][nlon]) v_field;
  float (*reference)[nclv][nlev][nlon] = (float (*)[nclv][nlev][nlon]) v_ref;

  zminval = +DBL_MAX;
  zmaxval = -DBL_MAX;
  zmaxerr = 0.0;
  zerrsum = 0.0;
  zsum = 0.0;
  zdiffmax = 0.0;

  // # pragma omp parallel for default(shared) private(b, bsize, jl, jk, jm) \
  //   reduction(min:zminval) reduction(max:zmaxval,zmaxerr) reduction(+:zerrsum,zsum)
  for (b = 0; b < nblocks; b++) {
    bsize = min(nlon, ngptot - b*nlon);  // field block size
    for (jm = 0; jm < nclv; jm++) {
      for (jl = 0; jl < nlev; jl++) {
	for (jk = 0; jk < bsize; jk++) {
          field_val = (double) field[b][jm][jl][jk];
          ref_val = (double) reference[b][jm][jl][jk];

	  zminval = fmin(zminval, field_val);
	  zmaxval = fmax(zmaxval, field_val);

	  // Difference against reference result in one-norm sense
	  zdiff = fabs(field_val - ref_val);
	  zmaxerr = fmax(zmaxerr, zdiff);
	  zdiffmax = fmax(zdiffmax, zdiff);
	  zerrsum = zerrsum + zdiff;
	  zsum = zsum + fabs(ref_val);
	}
      }
    }
  }
  zavgpgp = zerrsum / (double) ngptot;
  print_error(name, zminval, zmaxval, zmaxerr, zerrsum, zsum, zavgpgp, 3);
}


int cloudsc_validate(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma,
		     float *plude, float *pcovptot, float *prainfrac_toprfz, float *pfsqlf, float *pfsqif,
		     float *pfcqlng, float *pfcqnng, float *pfsqrf, float *pfsqsf, float *pfcqrng, float *pfcqsng,
		     float *pfsqltur, float *pfsqitur, float *pfplsl, float *pfplsn, float *pfhpsl, float *pfhpsn,
		     float *tend_loc_a, float *tend_loc_q, float *tend_loc_t, float *tend_loc_cld)
{
  const int nblocks = (ngptot / nproma) + min(ngptot % nproma, 1);
  float *ref_plude = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  float *ref_pcovptot = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  float *ref_prainfrac_toprfz = (float*) malloc( sizeof(float) * nblocks*nproma );
  float *ref_pfsqlf = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  float *ref_pfsqif = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  float *ref_pfcqlng = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  float *ref_pfcqnng = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  float *ref_pfsqrf = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  float *ref_pfsqsf = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  float *ref_pfcqrng = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  float *ref_pfcqsng = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  float *ref_pfsqltur = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  float *ref_pfsqitur = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  float *ref_pfplsl = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  float *ref_pfplsn = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  float *ref_pfhpsl = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  float *ref_pfhpsn = (float*) malloc( sizeof(float) * nblocks*(nlev+1)*nproma );
  float *ref_tend_loc_a = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  float *ref_tend_loc_q = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  float *ref_tend_loc_t = (float*) malloc( sizeof(float) * nblocks*nlev*nproma );
  float *ref_tend_loc_cld = (float*) malloc( sizeof(float) * nblocks*nclv*nlev*nproma );

  load_reference(nlon, nlev, nclv, ngptot, nproma,
                 ref_plude, ref_pcovptot, ref_prainfrac_toprfz, ref_pfsqlf, ref_pfsqif,
		 ref_pfcqlng, ref_pfcqnng, ref_pfsqrf, ref_pfsqsf, ref_pfcqrng, ref_pfcqsng,
		 ref_pfsqltur, ref_pfsqitur, ref_pfplsl, ref_pfplsn, ref_pfhpsl, ref_pfhpsn,
		 ref_tend_loc_a, ref_tend_loc_q, ref_tend_loc_t, ref_tend_loc_cld);


  printf(" %+20s %s %+20s %+20s %+20s %+20s %+20s\n",
	 "Variable", "Dim", "MinValue", "MaxValue", "AbsMaxErr", "AvgAbsErr/GP", "MaxRelErr-%");

  validate_2d("PLUDE", ref_plude, plude, nproma, nlev, ngptot, nblocks);
  validate_2d("PCOVPTOT", ref_pcovptot, pcovptot, nproma, nlev, ngptot, nblocks);
  validate_1d("PRAINFRAC_TOPRFZ", ref_prainfrac_toprfz, prainfrac_toprfz, nproma, ngptot, nblocks);
  validate_2d("PFSQLF", ref_pfsqlf, pfsqlf, nproma, nlev+1, ngptot, nblocks);
  validate_2d("PFSQIF", ref_pfsqif, pfsqif, nproma, nlev+1, ngptot, nblocks);
  validate_2d("PFCQLNG", ref_pfcqlng, pfcqlng, nproma, nlev+1, ngptot, nblocks);
  validate_2d("PFCQNNG", ref_pfcqnng, pfcqnng, nproma, nlev+1, ngptot, nblocks);
  validate_2d("PFSQRF", ref_pfsqrf, pfsqrf, nproma, nlev+1, ngptot, nblocks);
  validate_2d("PFSQSF", ref_pfsqsf, pfsqsf, nproma, nlev+1, ngptot, nblocks);
  validate_2d("PFCQRNG", ref_pfcqrng, pfcqrng, nproma, nlev+1, ngptot, nblocks);
  validate_2d("PFCQSNG", ref_pfcqsng, pfcqsng, nproma, nlev+1, ngptot, nblocks);
  validate_2d("PFSQLTUR", ref_pfsqltur, pfsqltur, nproma, nlev+1, ngptot, nblocks);
  validate_2d("PFSQITUR", ref_pfsqitur, pfsqitur, nproma, nlev+1, ngptot, nblocks);
  validate_2d("PFPLSL", ref_pfplsl, pfplsl, nproma, nlev+1, ngptot, nblocks);
  validate_2d("PFPLSN", ref_pfplsn, pfplsn, nproma, nlev+1, ngptot, nblocks);
  validate_2d("PFHPSL", ref_pfhpsl, pfhpsl, nproma, nlev+1, ngptot, nblocks);
  validate_2d("PFHPSN", ref_pfhpsn, pfhpsn, nproma, nlev+1, ngptot, nblocks);
  validate_2d("TENDENCY_LOC%A", ref_tend_loc_a, tend_loc_a, nproma, nlev, ngptot, nblocks);
  validate_2d("TENDENCY_LOC%Q", ref_tend_loc_q, tend_loc_q, nproma, nlev, ngptot, nblocks);
  validate_2d("TENDENCY_LOC%T", ref_tend_loc_t, tend_loc_t, nproma, nlev, ngptot, nblocks);
  validate_3d("TENDENCY_LOC%CLD", ref_tend_loc_cld, tend_loc_cld, nproma, nlev, nclv, ngptot, nblocks);

  free(ref_plude);
  free(ref_pcovptot);
  free(ref_prainfrac_toprfz);
  free(ref_pfsqlf);
  free(ref_pfsqif);
  free(ref_pfcqlng);
  free(ref_pfcqnng);
  free(ref_pfsqrf);
  free(ref_pfsqsf);
  free(ref_pfcqrng);
  free(ref_pfcqsng);
  free(ref_pfsqltur);
  free(ref_pfsqitur);
  free(ref_pfplsl);
  free(ref_pfplsn);
  free(ref_pfhpsl);
  free(ref_pfhpsn);
  free(ref_tend_loc_a);
  free(ref_tend_loc_q);
  free(ref_tend_loc_t);
  free(ref_tend_loc_cld);

  return 0;

}
