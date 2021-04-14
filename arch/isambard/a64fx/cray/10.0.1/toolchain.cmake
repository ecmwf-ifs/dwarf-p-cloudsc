####################################################################
# Compiler FLAGS
####################################################################

# General Flags (add to default)
set( ECBUILD_FIND_MPI ON )
set(CMAKE_POSITION_INDEPENDENT_CODE ON)

set(ECBUILD_Fortran_FLAGS_BIT "${ECBUILD_Fortran_FLAGS_BIT} -O3")
set(ECBUILD_Fortran_FLAGS_BIT "${ECBUILD_Fortran_FLAGS_BIT} -hfp1")
set(ECBUILD_Fortran_FLAGS_BIT "${ECBUILD_Fortran_FLAGS_BIT} -hscalar3")
set(ECBUILD_Fortran_FLAGS_BIT "${ECBUILD_Fortran_FLAGS_BIT} -hvector3")
set(ECBUILD_Fortran_FLAGS_BIT "${ECBUILD_Fortran_FLAGS_BIT} -G2")
set(ECBUILD_Fortran_FLAGS_BIT "${ECBUILD_Fortran_FLAGS_BIT} -hcontiguous")
set(ECBUILD_Fortran_FLAGS_BIT "${ECBUILD_Fortran_FLAGS_BIT} -haggress")
set(ECBUILD_Fortran_FLAGS_BIT "${ECBUILD_Fortran_FLAGS_BIT} -Ktrap=fp")

####################################################################
# LINK FLAGS
####################################################################

# Compatibility with HDF5 1.12
set(H5_USE_110_API ON)
