####################################################################
# COMPILER
####################################################################

include(CMakeForceCompiler)

set( ECBUILD_FIND_MPI OFF )

####################################################################
# COMMON FLAGS
####################################################################

set( ECBUILD_Fortran_FLAGS "\
-ffpe-trap=invalid,zero,overflow -fstack-arrays -fbacktrace -fno-second-underscore \
-ffree-form -ffast-math -fno-unsafe-math-optimizations  \
" )

set( ECBUILD_C_FLAGS "" )
