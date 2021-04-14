####################################################################
# Compiler FLAGS
####################################################################

# General Flags (add to default)
set( ECBUILD_FIND_MPI ON )
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -O3")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -hfp1")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -hscalar3")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -hvector3")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -G2")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -hcontiguous")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -haggress")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Ktrap=fp")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -DNDEBUG")

####################################################################
# LINK FLAGS
####################################################################

# Compatibility with HDF5 1.12
set(H5_USE_110_API ON)
