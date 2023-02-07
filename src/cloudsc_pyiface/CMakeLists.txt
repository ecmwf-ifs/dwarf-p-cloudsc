# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_FORTRAN_PYIFACE
    DESCRIPTION "Build the Python interface to the (slightly customized) cloudsc-fortran" DEFAULT ON
    CONDITION HDF5_FOUND
)
ecbuild_add_option( FEATURE CLOUDSC_FORTRAN_PYIFACE_BINARY
    DESCRIPTION "Build the debug-oriented binary for Python interface of (customized) cloudsc-fortran" DEFAULT ON 
    CONDITION HAVE_CLOUDSC_FORTRAN_PYIFACE 
)

if( HAVE_CLOUDSC_FORTRAN_PYIFACE )
   # Define module directory to facilitate f90wrap/f2py execution
   ecbuild_enable_fortran(MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../module)
   message (STATUS "Module directory set to:")
   message (STATUS  ${CMAKE_CURRENT_BINARY_DIR}/../../module)
   

   # Set names and lists to abstract the f90wrap/f2py process  
   set (DWARF_CLOUDSC_LIB  dwarf-cloudsc-lib)
   set (PYTHON_MODN cloudsc)
   set(FORTRAN_PYTHON_COMMON_SOURCES
            yomphyder
            yoecldp
            yoethf
            yomcst
            yoephli
            expand_mod
       )
   set(FORTRAN_PYTHON_DWARF_SOURCES
            cloudsc_driver_mod
       )
   set (COMMON_MOD_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/../common/module/)
   set (FORTRAN_SRC_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/fortransrc/)
   # Manipulate lists to prepare abstracted command-line input/output for f90wrap/f2py 
     set (F90WRAP_COMMON_SOURCES ${FORTRAN_PYTHON_COMMON_SOURCES})
     set ( F90WRAP_DWARF_SOURCES  ${FORTRAN_PYTHON_DWARF_SOURCES})
     set (   F2PY_COMMON_SOURCES ${FORTRAN_PYTHON_COMMON_SOURCES})
     set (    F2PY_DWARF_SOURCES  ${FORTRAN_PYTHON_DWARF_SOURCES})
     list(TRANSFORM F90WRAP_COMMON_SOURCES PREPEND ${COMMON_MOD_LOCATION})
     list(TRANSFORM F90WRAP_COMMON_SOURCES APPEND .F90 )
     list(TRANSFORM F90WRAP_DWARF_SOURCES  PREPEND ${FORTRAN_SRC_LOCATION})
     list(TRANSFORM F90WRAP_DWARF_SOURCES  APPEND .F90 )
     list(TRANSFORM F2PY_COMMON_SOURCES    PREPEND f90wrap_)
     list(TRANSFORM F2PY_COMMON_SOURCES    APPEND .f90 )
     list(TRANSFORM F2PY_DWARF_SOURCES     PREPEND f90wrap_)
     list(TRANSFORM F2PY_DWARF_SOURCES     APPEND .f90 )

   #Build CLOUDSC driver/kernel library, to be further linked by f2py
   ecbuild_add_library( TARGET ${DWARF_CLOUDSC_LIB}
        SOURCES
            ./fortransrc/cloudsc_driver_mod.F90
            ./fortransrc/cloudsc.F90
        LIBS
            cloudsc-common-lib
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
   )



if( HAVE_CLOUDSC_FORTRAN_PYIFACE_BINARY )
       # Define the (optional) binary build target for this variant
       ecbuild_add_executable( TARGET dwarf-cloudsc-fortran-pyiref
           SOURCES
               ./fortransrc/dwarf_cloudsc.F90
           LIBS
               cloudsc-common-lib   
                dwarf-cloudsc-lib
           DEFINITIONS ${CLOUDSC_DEFINITIONS}
       )


    ecbuild_add_test(
        TARGET dwarf-cloudsc-fortran-pyiref
        COMMAND bin/dwarf-cloudsc-fortran-pyiref
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
                   )
    endif()

    if( HAVE_HDF5 )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../input.h5 )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../reference.h5 )
    endif()

# F90wrap execution section

      add_custom_command(TARGET ${DWARF_CLOUDSC_LIB}  POST_BUILD
         COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
         COMMAND echo "Executing f90wrap from location: "
         COMMAND which f90wrap
         COMMAND  f90wrap -m${PYTHON_MODN} 
	                    ${F90WRAP_COMMON_SOURCES} 
			    ${F90WRAP_DWARF_SOURCES} 
			    -k ${CMAKE_CURRENT_SOURCE_DIR}/src/cloudscpytools/kind_map 
	   	 > f90wrap_log.txt 2> f90wrap_log_err.txt 
         COMMAND echo "Executing f2py-f90wrap from location: "
         COMMAND which f2py-f90wrap
         COMMAND ${CMAKE_COMMAND} -E env --unset=LD_LIBRARY_FLAGS 
	          LDFLAGS=-Wl,-rpath,../../lib NPY_DISTUTILS_APPEND_FLAGS=1  
		  f2py-f90wrap -c -m _${PYTHON_MODN}
	         -I${CMAKE_CURRENT_BINARY_DIR}/../../module
                 -I${CMAKE_CURRENT_BINARY_DIR}/../common/module
                 -L${CMAKE_CURRENT_BINARY_DIR}/../../../lib
                 -lcloudsc-common-lib 
		 -l${DWARF_CLOUDSC_LIB}
		   ${F2PY_COMMON_SOURCES} 
		   ${F2PY_DWARF_SOURCES} 
	     			     > f2py_log.txt 2> f2py_log_err.txt
          COMMAND echo "Runnning pip using the following Python executable"
          COMMAND echo ${Python_EXECUTABLE}
          COMMAND ${CMAKE_COMMAND} -E copy ${PYTHON_MODN}.py ./src/cloudscpytools
          COMMAND ${Python_EXECUTABLE} -m pip install -U .
          COMMAND ${CMAKE_COMMAND} -E make_directory ../../../bin/pythonexec
          COMMAND ${CMAKE_COMMAND} -E copy ${CMAKE_CURRENT_BINARY_DIR}/drivers/cloudsc_pydriver.py  ../../../bin/pythonexec/
          COMMENT "Completed execution of F90wrap/F2PY interface to cloudsc. "
          DEPENDS dwarf-cloudsc-lib
          VERBATIM
      )

          ecbuild_add_test(
             TARGET dwarf-cloudsc-lib
             COMMAND cloudsc_pydriver.py
             ARGS 1 100 16
             WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../../bin/pythonexec
             OMP 1
                          )
endif()
