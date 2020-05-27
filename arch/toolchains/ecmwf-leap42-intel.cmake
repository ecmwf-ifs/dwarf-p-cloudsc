####################################################################
# COMPILER
####################################################################

set( ECBUILD_FIND_MPI OFF )

####################################################################
# Compiler FLAGS
####################################################################

# General Flags (add to default)

set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -qopenmp-threadprivate compat")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -assume byterecl")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -convert big_endian")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -traceback")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -align array64byte")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -warn nounused,nouncalled")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -xHost")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -finline-functions")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -finline-limit=500")
set(ECBUILD_Fortran_FLAGS "${ECBUILD_Fortran_FLAGS} -Winline")
