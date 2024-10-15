####################################################################
# HIP 
####################################################################

set(CMAKE_HIP_ARCHITECTURES gfx90a)

####################################################################
# OpenMP
####################################################################

set( OpenMP_C_FLAGS           "-fopenmp" )
set( OpenMP_CXX_FLAGS         "-fopenmp" )
set( OpenMP_Fortran_FLAGS     "-fopenmp" )
set( OpenMP_C_LIB_NAMES       "craymp" )
set( OpenMP_CXX_LIB_NAMES     "craymp" )
set( OpenMP_Fortran_LIB_NAMES "craymp" )
set( OpenMP_craymp_LIBRARY    "craymp" )

####################################################################
# OpenACC
####################################################################

set( OpenACC_C_FLAGS       "-hacc" )
set( OpenACC_CXX_FLAGS     "-hacc" )
#set( OpenACC_Fortran_FLAGS "-hacc -hacc_model=auto_async_kernel:no_fast_addr:deep_copy" )
set( OpenACC_Fortran_FLAGS "-hacc" )

####################################################################
# General Flags (add to default)
####################################################################

set( ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -M2260" )

###################################################################
# Libraries
###################################################################

set( BLAS_LIBRARIES   "$ENV{CRAY_LIBSCI}" CACHE STRING "BLAS_LIBRARIES" FORCE )
set( LAPACK_LIBRARIES "$ENV{CRAY_LIBSCI}" CACHE STRING "LAPACK_LIBRARIES" FORCE )

#set ( CMAKE_C_FLAGS "${CMAKE_C_FLAGS} -craype-verbose" )
#set ( CMAKE_Fortran_FLAGS "${CMAKE_Fortran_FLAGS} -craype-verbose -fopenmp" )
#set ( CMAKE_EXE_LINKER_FLAGS  "${CMAKE_EXE_LINKER_FLAGS} -Wl, --as-needed" )
#set(CMAKE_SHARED_LINKER_FLAGS "${CMAKE_SHARED_LINKER_FLAGS} -fuse-ld=bfd")
#-set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-relax -lhugetlbfs")
#set(CMAKE_EXE_LINKER_FLAGS "${CMAKE_EXE_LINKER_FLAGS} -Wl,--no-relax -lhugetlbfs -fuse-ld=bfd -craype-verbose")
