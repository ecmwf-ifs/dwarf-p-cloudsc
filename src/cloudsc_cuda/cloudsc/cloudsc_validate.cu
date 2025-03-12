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

#include <math.h>
#include <float.h>

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

  printf(" %20s %dD%d %20.13le %20.13le %20.13le %20.13le %20.13le %s\n",
         name, ndim, iopt, zminval, zmaxval, zmaxerr, zavgpgp, zrelerr, clwarn);
}


void validate_1d(const char *name, dtype * v_ref, dtype * v_field, int nlon, int ngptot, int nblocks)
{
  /* Computes and prints errors in the "L2 norm sense" */
  int b, bsize, jk;
  double zminval, zmaxval, zdiff, zmaxerr, zerrsum, zsum, zrelerr, zavgpgp;
  double field_val, ref_val;
  dtype (*field)[nlon] = (dtype (*)[nlon]) v_field;
  dtype (*reference)[nlon] = (dtype (*)[nlon]) v_ref;

  zminval = +DBL_MAX;
  zmaxval = -DBL_MAX;
  zmaxerr = 0.0;
  zerrsum = 0.0;
  zsum = 0.0;

  # pragma omp parallel for default(shared) private(b, bsize, jk)		\
     reduction(min:zminval) reduction(max:zmaxval,zmaxerr) reduction(+:zerrsum,zsum)
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


void validate_2d(const char *name, dtype *v_ref, dtype *v_field, int nlon, int nlev, int ngptot, int nblocks)
{
  /* Computes and prints errors in the "L2 norm sense" */
  int b, bsize, jl, jk;
  double zminval, zmaxval, zdiff, zmaxerr, zerrsum, zsum, zrelerr, zavgpgp, zdiffmax;
  double field_val, ref_val;
  dtype (*field)[nlev][nlon] = (dtype (*)[nlev][nlon]) v_field;
  dtype (*reference)[nlev][nlon] = (dtype (*)[nlev][nlon]) v_ref;

  zminval = +DBL_MAX;
  zmaxval = -DBL_MAX;
  zmaxerr = 0.0;
  zerrsum = 0.0;
  zsum = 0.0;
  zdiffmax = 0.0;

  # pragma omp parallel for default(shared) private(b, bsize, jl, jk) \
     reduction(min:zminval) reduction(max:zmaxval,zmaxerr) reduction(+:zerrsum,zsum)
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


void validate_3d(const char *name, dtype *v_ref, dtype *v_field, int nlon,
    int nlev, int nclv, int ngptot, int nblocks)
{
  /* Computes and prints errors in the "L2 norm sense" */
  int b, bsize, jl, jk, jm;
  double zminval, zmaxval, zdiff, zmaxerr, zerrsum, zsum, zrelerr, zavgpgp, zdiffmax;
  double field_val, ref_val;
  dtype (*field)[nclv][nlev][nlon] = (dtype (*)[nclv][nlev][nlon]) v_field;
  dtype (*reference)[nclv][nlev][nlon] = (dtype (*)[nclv][nlev][nlon]) v_ref;

  zminval = +DBL_MAX;
  zmaxval = -DBL_MAX;
  zmaxerr = 0.0;
  zerrsum = 0.0;
  zsum = 0.0;
  zdiffmax = 0.0;

  #pragma omp parallel for default(shared) private(b, bsize, jl, jk, jm) \
     reduction(min:zminval) reduction(max:zmaxval,zmaxerr) reduction(+:zerrsum,zsum)
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
		     dtype *plude, dtype *pcovptot, dtype *prainfrac_toprfz, dtype *pfsqlf, dtype *pfsqif,
		     dtype *pfcqlng, dtype *pfcqnng, dtype *pfsqrf, dtype *pfsqsf, dtype *pfcqrng, dtype *pfcqsng,
		     dtype *pfsqltur, dtype *pfsqitur, dtype *pfplsl, dtype *pfplsn, dtype *pfhpsl, dtype *pfhpsn,
		     dtype *tend_loc_a, dtype *tend_loc_q, dtype *tend_loc_t, dtype *tend_loc_cld)
{
  const int nblocks = (ngptot / nproma) + min(ngptot % nproma, 1);
  dtype *ref_plude = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  dtype *ref_pcovptot = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  dtype *ref_prainfrac_toprfz = (dtype*) malloc( sizeof(dtype) * nblocks*nproma );
  dtype *ref_pfsqlf = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  dtype *ref_pfsqif = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  dtype *ref_pfcqlng = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  dtype *ref_pfcqnng = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  dtype *ref_pfsqrf = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  dtype *ref_pfsqsf = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  dtype *ref_pfcqrng = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  dtype *ref_pfcqsng = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  dtype *ref_pfsqltur = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  dtype *ref_pfsqitur = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  dtype *ref_pfplsl = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  dtype *ref_pfplsn = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  dtype *ref_pfhpsl = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  dtype *ref_pfhpsn = (dtype*) malloc( sizeof(dtype) * nblocks*(nlev+1)*nproma );
  dtype *ref_tend_loc_a = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  dtype *ref_tend_loc_q = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  dtype *ref_tend_loc_t = (dtype*) malloc( sizeof(dtype) * nblocks*nlev*nproma );
  dtype *ref_tend_loc_cld = (dtype*) malloc( sizeof(dtype) * nblocks*nclv*nlev*nproma );

  load_reference(nlon, nlev, nclv, ngptot, nproma,
                 ref_plude, ref_pcovptot, ref_prainfrac_toprfz, ref_pfsqlf, ref_pfsqif,
		 ref_pfcqlng, ref_pfcqnng, ref_pfsqrf, ref_pfsqsf, ref_pfcqrng, ref_pfcqsng,
		 ref_pfsqltur, ref_pfsqitur, ref_pfplsl, ref_pfplsn, ref_pfhpsl, ref_pfhpsn,
		 ref_tend_loc_a, ref_tend_loc_q, ref_tend_loc_t, ref_tend_loc_cld);


  printf(" %20s %s %20s %20s %20s %20s %20s\n",
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
