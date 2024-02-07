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
#include <limits>
//#include <limits.h>


#define min(a,b) (((a)<(b))?(a):(b))

void print_error(const char *name, double zminval, double zmaxval, double zmaxerr,
		 double zerrsum, double zsum, double zavgpgp, int ndim)
{
  double zrelerr, zeps = std::numeric_limits<double>::epsilon();
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

  printf(" %+20s %dD%d %20.13le %20.13le %20.13le %20.13le %20.13le %s\n",
         name, ndim, iopt, zminval, zmaxval, zmaxerr, zavgpgp, zrelerr, clwarn);
}


void validate_1d(const char *name, double * v_ref, double * v_field, int nlon, int ngptot, int nblocks)
{
  /* Computes and prints errors in the "L2 norm sense" */
  int b, bsize, jk;
  double zminval, zmaxval, zdiff, zmaxerr, zerrsum, zsum, zrelerr, zavgpgp;
  //double (*field)[nlon] = (double (*)[nlon]) v_field;
  //double (*reference)[nlon] = (double (*)[nlon]) v_ref;

  zminval = +std::numeric_limits<double>::max();
  zmaxval = -std::numeric_limits<double>::max();
  zmaxerr = 0.0;
  zerrsum = 0.0;
  zsum = 0.0;

//  #pragma omp parallel for default(shared) private(b, bsize, jk)		\
//    reduction(min:zminval) reduction(max:zmaxval,zmaxerr) reduction(+:zerrsum,zsum)
  for (b = 0; b < nblocks; b++) {
    bsize = min(nlon, ngptot - b*nlon);  // field block size
    for (jk = 0; jk < bsize; jk++) {
      zminval = fmin(zminval, v_field[b*nlon+jk]);
      zmaxval = fmax(zmaxval, v_field[b*nlon+jk]);

      // Difference against reference result in one-norm sense
      zdiff = fabs(v_field[b*nlon+jk] - v_ref[b*nlon+jk]);
      zmaxerr = fmax(zmaxerr, zdiff);
      zerrsum = zerrsum + zdiff;
      zsum = zsum + abs(v_ref[b*nlon+jk]);
    }
  }
  zavgpgp = zerrsum / (double) ngptot;
  print_error(name, zminval, zmaxval, zmaxerr, zerrsum, zsum, zavgpgp, 2);
}


void validate_2d(const char *name, double *v_ref, double *v_field, int nlon, int nlev, int ngptot, int nblocks)
{
  /* Computes and prints errors in the "L2 norm sense" */
  int b, bsize, jl, jk;
  double zminval, zmaxval, zdiff, zmaxerr, zerrsum, zsum, zrelerr, zavgpgp;
//  double (*field)[nlev][nlon] = (double (*)[nlev][nlon]) v_field;
//  double (*reference)[nlev][nlon] = (double (*)[nlev][nlon]) v_ref;

  zminval = +std::numeric_limits<double>::max();
  zmaxval = -std::numeric_limits<double>::max();
  zmaxerr = 0.0;
  zerrsum = 0.0;
  zsum = 0.0;

//  #pragma omp parallel for default(shared) private(b, bsize, jl, jk) \
//    reduction(min:zminval) reduction(max:zmaxval,zmaxerr) reduction(+:zerrsum,zsum)

  for (b = 0; b < nblocks; b++) {
    bsize = min(nlon, ngptot - b*nlon);  // field block size
    for (jl = 0; jl < nlev; jl++) {
      for (jk = 0; jk < bsize; jk++) {
	zminval = fmin(zminval, v_field[b*nlev*nlon+jl*nlon+jk]);
	zmaxval = fmax(zmaxval, v_field[b*nlev*nlon+jl*nlon+jk]);

	// Difference against reference result in one-norm sense
	zdiff = fabs(v_field[b*nlev*nlon+jl*nlon+jk] - v_ref[b*nlev*nlon+jl*nlon+jk]);
	zmaxerr = fmax(zmaxerr, zdiff);
	zerrsum = zerrsum + zdiff;
	zsum = zsum + abs(v_ref[b*nlev*nlon+jl*nlon+jk]);
      }
    }
  }

  zavgpgp = zerrsum / (double) ngptot;
  print_error(name, zminval, zmaxval, zmaxerr, zerrsum, zsum, zavgpgp, 2);
}


void validate_3d(const char *name, double *v_ref, double *v_field, int nlon,
    int nlev, int nclv, int ngptot, int nblocks)
{
  /* Computes and prints errors in the "L2 norm sense" */
  int b, bsize, jl, jk, jm;
  double zminval, zmaxval, zdiff, zmaxerr, zerrsum, zsum, zrelerr, zavgpgp;
//  double (*field)[nclv][nlev][nlon] = (double (*)[nclv][nlev][nlon]) v_field;
//  double (*reference)[nclv][nlev][nlon] = (double (*)[nclv][nlev][nlon]) v_ref;

  zminval = +std::numeric_limits<double>::max();
  zmaxval = -std::numeric_limits<double>::max();
  zmaxerr = 0.0;
  zerrsum = 0.0;
  zsum = 0.0;

//  #pragma omp parallel for default(shared) private(b, bsize, jl, jk, jm) \
//    reduction(min:zminval) reduction(max:zmaxval,zmaxerr) reduction(+:zerrsum,zsum)
  for (b = 0; b < nblocks; b++) {
    bsize = min(nlon, ngptot - b*nlon);  // field block size
    for (jm = 0; jm < nclv; jm++) {
      for (jl = 0; jl < nlev; jl++) {
	for (jk = 0; jk < bsize; jk++) {
	  zminval = fmin(zminval, v_field[b*nclv*nlev*nlon+jm*nlev*nlon+jl*nlon+jk]);
	  zmaxval = fmax(zmaxval, v_field[b*nclv*nlev*nlon+jm*nlev*nlon+jl*nlon+jk]);

	  // Difference against reference result in one-norm sense
	  zdiff = fabs(v_field[b*nclv*nlev*nlon+jm*nlev*nlon+jl*nlon+jk] - v_ref[b*nclv*nlev*nlon+jm*nlev*nlon+jl*nlon+jk]);
	  zmaxerr = fmax(zmaxerr, zdiff);
	  zerrsum = zerrsum + zdiff;
	  zsum = zsum + abs(v_ref[b*nclv*nlev*nlon+jm*nlev*nlon+jl*nlon+jk]);
	}
      }
    }
  }
  zavgpgp = zerrsum / (double) ngptot;
  print_error(name, zminval, zmaxval, zmaxerr, zerrsum, zsum, zavgpgp, 2);
}


int cloudsc_validate(const int nlon, const int nlev, const int nclv, const int ngptot, const int nproma,
		     double *plude, double *pcovptot, double *prainfrac_toprfz, double *pfsqlf, double *pfsqif,
		     double *pfcqlng, double *pfcqnng, double *pfsqrf, double *pfsqsf, double *pfcqrng, double *pfcqsng,
		     double *pfsqltur, double *pfsqitur, double *pfplsl, double *pfplsn, double *pfhpsl, double *pfhpsn,
		     double *tend_loc_a, double *tend_loc_q, double *tend_loc_t, double *tend_loc_cld)
{
  const int nblocks = (ngptot / nproma) + min(ngptot % nproma, 1);
  double *ref_plude = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  double *ref_pcovptot = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  double *ref_prainfrac_toprfz = (double*) malloc( sizeof(double) * nblocks*nproma );
  double *ref_pfsqlf = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  double *ref_pfsqif = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  double *ref_pfcqlng = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  double *ref_pfcqnng = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  double *ref_pfsqrf = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  double *ref_pfsqsf = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  double *ref_pfcqrng = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  double *ref_pfcqsng = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  double *ref_pfsqltur = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  double *ref_pfsqitur = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  double *ref_pfplsl = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  double *ref_pfplsn = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  double *ref_pfhpsl = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  double *ref_pfhpsn = (double*) malloc( sizeof(double) * nblocks*(nlev+1)*nproma );
  double *ref_tend_loc_a = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  double *ref_tend_loc_q = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  double *ref_tend_loc_t = (double*) malloc( sizeof(double) * nblocks*nlev*nproma );
  double *ref_tend_loc_cld = (double*) malloc( sizeof(double) * nblocks*nclv*nlev*nproma );

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
