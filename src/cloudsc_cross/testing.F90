! (C) Copyright 1988- ECMWF.
!
! This software is licensed under the terms of the Apache Licence Version 2.0
! which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
!
! In applying this licence, ECMWF does not waive the privileges and immunities
! granted to it by virtue of its status as an intergovernmental organisation
! nor does it submit to any jurisdiction.

PROGRAM DWARF_CLOUDSC


#if GPU_OFFLOAD == ACC_OFFLOAD
print *, "offload via ACC"
#elif GPU_OFFLOAD == OMP_OFFLOAD
print *, "offload via OMP"
#else
print *, "offload undefined"
#endif

! print *, "GPU_LANG: ", GPU_LANG
! print *, "CUDA LANG: ", CUDA_LANG
! print *, "HIP_LANG: ", HIP_LANG
! print *, "SYCL LANG: ", SYCL_LANG

#if GPU_LANG == CUDA_LANG
print *, "GPU lang CUDA"
#elif GPU_LANG == HIP_LANG
print *, "GPU lang HIP"
#elif GPU_LANG == SYCL_LANG
print *, "GPU lang SYCL"
#else
print * "GPU lang undefined"
#endif

! print *, "Hi"
! print *, "OFFLOAD: ", OFFLOAD
! print *, "CROSS_OFFLOAD: ", CROSS_OFFLOAD
! print *, "CLOUDSC_CROSS_OFFLOAD: ", CLOUDSC_CROSS_OFFLOAD
! print *, "CLOUDSC_CROSS_LANG:    ", CLOUDSC_CROSS_LANG

! #ifdef CLOUDSC_CROSS_OFFLOAD
! print *, "offload is ACC"
! #else
! print *, "offload not ACC"
! #endif

END PROGRAM DWARF_CLOUDSC
