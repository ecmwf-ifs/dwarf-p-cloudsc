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
    cloudsc_driver(omp_threads, ngptot, nproma);
  }
  else {
    printf("Calling c-cloudsc with the right number of arguments will work better ;-) \n",argc);
    return_code = EXIT_FAILURE;
  }

   return return_code;
}
