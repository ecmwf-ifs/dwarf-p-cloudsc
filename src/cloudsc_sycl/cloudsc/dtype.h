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

// #define MYMAX(x,y) fmaxf(x,y)
// #define MYMIN(x,y) fminf(x,y)
// #define MYEXP(x) exp(x)
// #define MYPOW(x,y) pow(x,y)
// #define MYPOWN(x,y) pow(x,y)
// #define MYABS(x) fabs(x)

#if GPU_MATH == C_MATH

#define MYMAX(x,y) fmax(x,y)
#define MYMIN(x,y) fmin(x,y)
#define MYEXP(x) exp(x)
#define MYPOW(x,y) pow(x,y)
#define MYPOWN(x,y) pow(x,(float)y)
#define MYABS(x) fabs(x)

#elif GPU_MATH == STD_MATH

#define MYMAX(x,y) std::max((float)x,(float)y)
#define MYMIN(x,y) std::min((float)x,(float)y)
#define MYEXP(x) std::exp(x)
#define MYPOW(x,y) std::pow(x,y)
#define MYPOWN(x,y) std::pow(x,y)
#define MYABS(x) std::abs(x)

#elif GPU_MATH == SYCL_MATH

#define MYMAX(x,y) cl::sycl::max(x,y)
#define MYMIN(x,y) cl::sycl::min(x,y)
#define MYEXP(x) cl::sycl::exp(x)
#define MYPOW(x,y) cl::sycl::pow(x,y)
#define MYPOWN(x,y) cl::sycl::pown(x,y)
#define MYABS(x) cl::sycl::abs(x)

#endif // GPU_MATH

#else

// #define MYMAX(x,y) fmax(x,y)
// #define MYMIN(x,y) fmin(x,y)
// #define MYEXP(x) exp(x)
// #define MYPOW(x,y) pow(x,y)
// #define MYPOWN(x,y) pow(x,y)
// #define MYABS(x) fabs(x)

#if GPU_MATH == C_MATH

#define MYMAX(x,y) fmax(x,y)
#define MYMIN(x,y) fmin(x,y)
#define MYEXP(x) exp(x)
#define MYPOW(x,y) pow(x,y)
#define MYPOWN(x,y) pow(x,y)
#define MYABS(x) fabs(x)

#elif GPU_MATH == STD_MATH

#define MYMAX(x,y) std::max(x,y)
#define MYMIN(x,y) std::min(x,y)
#define MYEXP(x) std::exp(x)
#define MYPOW(x,y) std::pow(x,y)
#define MYPOWN(x,y) std::pow(x,y)
#define MYABS(x) std::abs(x)

#elif GPU_MATH == SYCL_MATH

#define MYMAX(x,y) cl::sycl::max(x,y)
#define MYMIN(x,y) cl::sycl::min(x,y)
#define MYEXP(x) cl::sycl::exp(x)
#define MYPOW(x,y) cl::sycl::pow(x,y)
#define MYPOWN(x,y) cl::sycl::pown(x,y)
#define MYABS(x) cl::sycl::abs(x)

#endif // GPU_MATH
#endif // SINGLE

#endif // CLOUDSC_DTYPE_H
