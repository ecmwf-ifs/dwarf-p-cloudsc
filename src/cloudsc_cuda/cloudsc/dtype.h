/*
 * (C) Copyright 1988- ECMWF.
 *
 * This software is licensed under the terms of the Apache Licence Version 2.0
 * which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
 * In applying this licence, ECMWF does not waive the privileges and immunities
 * granted to it by virtue of its status as an intergovernmental organisation
 * nor does it submit to any jurisdiction.
 */

#ifndef CLOUDSC_DTYPE_H
#define CLOUDSC_DTYPE_H

// FLOATING POINT PRECISION
#ifdef SINGLE
typedef float dtype;
#else
typedef double dtype;
#endif


// MATH FUNCTIONS
#ifdef SINGLE

#define MYMAX(x,y) fmaxf(x,y)
#define MYMIN(x,y) fminf(x,y)
#define MYEXP(x) exp(x)
#define MYPOW(x,y) pow(x,y)
#define MYPOWN(x,y) pow(x,y)
#define MYABS(x) fabs(x)

#else

#define MYMAX(x,y) fmax(x,y)
#define MYMIN(x,y) fmin(x,y)
#define MYEXP(x) exp(x)
#define MYPOW(x,y) pow(x,y)
#define MYPOWN(x,y) pow(x,y)
#define MYABS(x) fabs(x)

#endif

#endif // CLOUDSC_DTYPE_H
