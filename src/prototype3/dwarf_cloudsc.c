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
