#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#include "load_state.h"
#include "cloudsc_c.h"

void cloudsc_driver(int numthreads, int numcols, int nproma);
