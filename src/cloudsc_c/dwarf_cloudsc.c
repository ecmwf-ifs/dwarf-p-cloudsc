/*
 * (C) Copyright 1988- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <stdio.h>
#include <stdlib.h>
#ifdef _OPENMP
#include <omp.h>
#endif

#include "cloudsc_driver.h"


int main( int argc, char *argv[] ) {

  int omp_threads, ngptot, nproma;
  int return_code;

  return_code = 0;

  // default values
  omp_threads = 1;
  ngptot = 100;
  nproma = 4;

  if (argc == 1) {
    cloudsc_driver(omp_threads, ngptot, nproma);
  }
  else if (argc == 4) {
    omp_threads = atoi( argv[1] );
    ngptot      = atoi( argv[2] );
    nproma      = atoi( argv[3] );
    if (omp_threads <= 0) {
#ifdef _OPENMP
      omp_threads = omp_get_max_threads();
#else
      // if arg is 0 or negative, and OpenMP disabled; defaults to 1
      omp_threads = 1;
#endif
    }
    cloudsc_driver(omp_threads, ngptot, nproma);
  }
  else {
    printf("Calling c-cloudsc with the right number of arguments will work better ;-) \n");
    return_code = EXIT_FAILURE;
  }

   return return_code;
}
