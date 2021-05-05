####################################################################
# COMPILER
####################################################################

set( ECBUILD_FIND_MPI ON )
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

if( ENABLE_MPI )
set(CMAKE_Fortran_COMPILER "mpifrt")
set(CMAKE_C_COMPILER "mpifcc")
set(CMAKE_CXX_COMPILER "mpiFCC")
endif()

####################################################################
# OpenMP FLAGS
####################################################################

set( OpenMP_C_FLAGS             "-Kopenmp" )
set( OpenMP_CXX_FLAGS           "-Kopenmp" )
set( OpenMP_Fortran_FLAGS       "-Kopenmp" )

####################################################################
# COMMON FLAGS
####################################################################

set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -fpic")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Kfast")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Kopenmp")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Nfjomplib")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -O3")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -KA64FX")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -KSVE")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -KARMV8_3_A")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Ksimd=2")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Klto")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -SSL2")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Kassume=notime_saving_compilation")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -DNDEBUG")

set(ECBUILD_C_FLAGS "${ECBUILD_C_FLAGS} -fpic")
set(ECBUILD_C_FLAGS "${ECBUILD_C_FLAGS} -Klto")
set(ECBUILD_CXX_FLAGS "${ECBUILD_CXX_FLAGS} -fpic")
set(ECBUILD_CXX_FLAGS "${ECBUILD_CXX_FLAGS} -Klto")

set(ECBUILD_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS} -Kopenmp")
set(ECBUILD_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS} -Nfjomplib")
set(ECBUILD_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS} -Klto")
set(ECBUILD_LINKER_FLAGS "${ECBUILD_LINKER_FLAGS} -SSL2")

# Compatibility with HDF5 1.12
set(H5_USE_110_API ON)
