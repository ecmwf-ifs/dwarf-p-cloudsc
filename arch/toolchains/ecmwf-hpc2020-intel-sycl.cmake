# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.



####################################################################
# Compiler FLAGS
####################################################################

# General Flags (add to default)

set(ECBUILD_Fortran_FLAGS "-g")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -qopenmp-threadprivate compat")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -assume byterecl")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -convert big_endian")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -traceback")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -align array64byte")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -warn nounused,nouncalled")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -march=core-avx2")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -finline-functions")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -finline-limit=1500")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Winline")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -no-fma")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -assume realloc_lhs")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fp-model precise")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -ftz")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fp-speculation=safe")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fast-transcendentals")

####################################################################
# Additional compiler flags for SYCL offload via CUDA backend
####################################################################

# Additional Intel DPCPP compiler for SYCL offload
set(CMAKE_CXX_COMPILER "/home/nams/opt/dpcpp/bin/clang++")

# Initial set of flags to things going with a custom DPCPP install on AC
set(CMAKE_CXX_FLAGS "-O3 -L/home/nams/opt/dpcpp/lib -fopenmp -lstdc++")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsycl-early-optimizations -fsycl")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -fsycl-targets=nvptx64-nvidia-cuda")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -Xsycl-target-backend --cuda-gpu-arch=sm_80")
set(CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS} -I/usr/local/apps/intel/2021.4.0/compiler/2021.4.0/linux/compiler/include")

