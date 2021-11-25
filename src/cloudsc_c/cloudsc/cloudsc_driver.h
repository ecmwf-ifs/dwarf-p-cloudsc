/*
 * (C) Copyright 1988- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#include <stdlib.h>
#include <stdio.h>
#include <string.h>
#include <omp.h>

#include "yomcst_c.h"
#include "yoethf_c.h"
#include "yoecldp_c.h"

#include "load_state.h"
#include "cloudsc_c.h"
#include "cloudsc_validate.h"

void cloudsc_driver(int numthreads, int numcols, int nproma);
