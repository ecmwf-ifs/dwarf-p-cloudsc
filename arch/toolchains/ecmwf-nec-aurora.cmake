
####################################################################
# COMPILER
####################################################################

include( /opt/nec/ve/share/cmake/toolchainVE.cmake )


set( ECBUILD_FIND_MPI ON )

####################################################################
# Enviroment Variables
####################################################################
set(NMPI_ROOT /opt/nec/ve/mpi/2.23.0)

####################################################################
# OpenMP FLAGS
####################################################################

set( OpenMP_C_FLAGS             "-fopenmp " )
set( OpenMP_CXX_FLAGS           "-fopenmp " )
set( OpenMP_Fortran_FLAGS       "-fopenmp " )

####################################################################
# OpenAcc FLAGS
####################################################################

set( OpenACC_Fortran_FLAGS "-acc -ta=tesla:lineinfo,deepcopy,maxregcount:100,fastmath" )
set( OpenACC_Fortran_FLAGS "${OpenACC_Fortran_FLAGS} -Mvect=levels:6" )
set( OpenACC_Fortran_FLAGS "${OpenACC_Fortran_FLAGS} -Mconcur=levels:6" )
set( OpenACC_Fortran_FLAGS "${OpenACC_Fortran_FLAGS} -Minfo" )

####################################################################
# NEC MPI Compiler
####################################################################

set(MPI_C_COMPILER ${NMPI_ROOT}/bin/mpincc CACHE FILEPATH "")
set(MPI_C_INCLUDE_PATH ${NMPI_ROOT}/include CACHE FILEPATH "")
set(MPI_C_LIBRARIES ${NMPI_ROOT}/lib64/ve/libmpi.a CACHE FILEPATH "")
set(MPI_C_COMPILE_FLAGS "-D_MPIPP_INCLUDE" CACHE STRING "")

set(MPI_CXX_COMPILER ${NMPI_ROOT}/bin/mpinc++ CACHE FILEPATH "")
set(MPI_CXX_INCLUDE_PATH ${NMPI_ROOT}/include CACHE FILEPATH "")
set(MPI_CXX_LIBRARIES ${NMPI_ROOT}/lib64/ve/libmpi++.a CACHE FILEPATH "")

set(MPI_Fortran_COMPILER ${NMPI_ROOT}/bin/mpifort CACHE FILEPATH "")
set(MPI_Fortran_INCLUDE_PATH ${NMPI_ROOT}/include CACHE FILEPATH "")
set(MPI_Fortran_ADDITIONAL_INCLUDE_DIR ${NMPI_ROOT}/lib/ve/module CACHE FILEPATH "")
set(MPI_Fortran_LIBRARIES ${NMPI_ROOT}/lib64/ve/libmpi.a CACHE FILEPATH "")
set(MPI_Fortran_COMPILE_FLAGS "-D_MPIPP_INCLUDE" CACHE STRING "")
####################################################################
# COMMON FLAGS
####################################################################

set(ECBUILD_Fortran_FLAGS "-fpic")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -mstack-arrays")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fdiag-vector=3")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fcse-after-vectorization")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -floop-collapse ")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -floop-fusion ")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -floop-interchange ")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -floop-unroll-complete=200 ")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -ftrace")
###set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fmove-loop-invariants-if ")
###set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -freplace-loop-equation ")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -msched-interblock ")
###set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -mvector-floating-divide-instruction ")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -mvector-power-to-explog ")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -mvector-sqrt-instruction ")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -mvector-threshold=3 ")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -finline-functions ")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -finline-max-depth=5 ")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -finline-max-function-size=200 ")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -mvector-merge-conditional ")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fivdep ")
###set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -floop-strip-mine ")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -muse-mmap ")
##set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -mvector-packed")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -report-all")

set( ECBUILD_Fortran_FLAGS_BIT "-O4 -mvector-fma" )

set( ECBUILD_C_FLAGS "-O2 " )

set( ECBUILD_CXX_FLAGS "-O2" )

# Fix for C++ template headers needed for Serialbox
set( GNU_HEADER_INCLUDE "-I/usr/local/apps/gcc/7.3.0/lib/gcc/x86_64-linux-gnu/7.3.0/include-fixed" )
set( ECBUILD_CXX_FLAGS "${ECBUILD_CXX_FLAGS} ${GNU_HEADER_INCLUDE}" )
