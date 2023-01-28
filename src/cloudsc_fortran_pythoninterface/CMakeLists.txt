# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_FORTRAN_PYIFACE
    DESCRIPTION "Build the Python interface to cloudsc-fortran" DEFAULT ON
    CONDITION HDF5_FOUND
)

if( HAVE_CLOUDSC_FORTRAN_PYIFACE )
   message (STATUS "Setting up module directory")
   ecbuild_enable_fortran(MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../module)
   message (STATUS  ${CMAKE_CURRENT_BINARY_DIR}/../../module)

    set (DWARF_CLOUDSC_LIB  dwarf-cloudsc-lib)
    ecbuild_add_library( TARGET ${DWARF_CLOUDSC_LIB}
        SOURCES
            ../cloudsc_fortran2/cloudsc_driver_mod.F90
            ../cloudsc_fortran2/cloudsc.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )
    set (PYTHON_MODN cloudsc)

     set(FORTRAN_PYTHON_COMMON_SOURCES
            yomphyder
            yoecldp
            yoethf
            yomcst
            yoephli
        )

     set(FORTRAN_PYTHON_DWARF_SOURCES
            cloudsc_driver_mod
        )

     set (F90WRAP_COMMON_SOURCES ${FORTRAN_PYTHON_COMMON_SOURCES})
     set (F90WRAP_DWARF_SOURCES  ${FORTRAN_PYTHON_DWARF_SOURCES})
     set (F2PY_COMMON_SOURCES ${FORTRAN_PYTHON_COMMON_SOURCES})
     set (F2PY_DWARF_SOURCES  ${FORTRAN_PYTHON_DWARF_SOURCES})
     list(TRANSFORM F90WRAP_COMMON_SOURCES PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/../common/module/)
     list(TRANSFORM F90WRAP_COMMON_SOURCES APPEND .F90 )
     list(TRANSFORM F90WRAP_DWARF_SOURCES PREPEND ${CMAKE_CURRENT_SOURCE_DIR}/../cloudsc_fortran2/)
     list(TRANSFORM F90WRAP_DWARF_SOURCES APPEND .F90 )
     list(TRANSFORM F2PY_COMMON_SOURCES PREPEND f90wrap_)
     list(TRANSFORM F2PY_COMMON_SOURCES APPEND .f90 )
     list(TRANSFORM F2PY_DWARF_SOURCES PREPEND f90wrap_)
     list(TRANSFORM F2PY_DWARF_SOURCES APPEND .f90 )

  set(PYTHON_SOURCES
            ./dwarfdriver.py
            ./cloudsc_data.py
        )

   add_custom_target( pythonexecs ALL SOURCES ${PYTHON_SOURCES}
#        COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}/pythonsrc  ${CMAKE_CURRENT_BINARY_DIR}/pythonsrc
         COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}  ../../bin/pythonexec
         COMMAND sh -c "rm -f ../../bin/pythonexec/*.f90"
         COMMENT "Copy Python driver to Fortran kernel to executable directory"
        )

    target_link_libraries( ${DWARF_CLOUDSC_LIB} PRIVATE cloudsc-common-lib HDF5::HDF5 )
    # Define the binary build target for this variant
    ecbuild_add_executable( TARGET dwarf-cloudsc-fortran-pyiref
        SOURCES
            ../cloudsc_fortran2/dwarf_cloudsc.F90
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    target_link_libraries( dwarf-cloudsc-fortran-pyiref PRIVATE cloudsc-common-lib ${DWARF_CLOUDSC_LIB})

    target_include_directories(dwarf-cloudsc-fortran-pyiref PRIVATE ${CMAKE_CURRENT_BINARY_DIR})
    target_include_directories(dwarf-cloudsc-fortran-pyiref PRIVATE ${CMAKE_CURRENT_BINARY_DIR}/../common)

    if( HAVE_OMP AND TARGET OpenMP::OpenMP_Fortran )
        target_link_libraries( dwarf-cloudsc-fortran-pyiref PRIVATE OpenMP::OpenMP_Fortran )
    endif()

    if( HAVE_MPI AND TARGET MPI::MPI_Fortran )
        target_link_libraries( dwarf-cloudsc-fortran-pyiref PRIVATE MPI::MPI_Fortran )
    endif()

    # Create symlink for the input data
    if( HAVE_SERIALBOX )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )
    endif()

    if( HAVE_HDF5 )
      target_include_directories( dwarf-cloudsc-fortran-pyiref PRIVATE ${HDF5_Fortran_INCLUDE_DIRS} )
      target_link_libraries( dwarf-cloudsc-fortran-pyiref PRIVATE ${HDF5_LIBRARIES} )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../input.h5 )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../reference.h5 )
    endif()

    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../input.h5 )
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../reference.h5 )

# F90wrap execution section
    execute_process(COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR}  ${CMAKE_CURRENT_BINARY_DIR})
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/input.h5 )
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
      ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/reference.h5 )

      add_custom_command(TARGET dwarf-cloudsc-fortran-pyiref POST_BUILD
         COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_SOURCE_DIR} ${CMAKE_CURRENT_BINARY_DIR}
	 #COMMAND python3 -m venv cstest
	 #COMMAND source ./cstest/bin/activate
         COMMAND pip install --upgrade pip
	 COMMAND pip install f90wrap h5py
         COMMAND echo "Executing f90wrap"
         COMMAND  f90wrap -m${PYTHON_MODN} ${F90WRAP_COMMON_SOURCES} ${F90WRAP_DWARF_SOURCES} -k ${CMAKE_CURRENT_SOURCE_DIR}/kind_map -P 
	 #> f90wrap_log.txt 2> f90wrap_log_err.txt COMMAND echo "Executing f2py-f90wrap"
         COMMAND pwd 
         COMMAND  f2py-f90wrap -c -m _${PYTHON_MODN} --include-paths ${CMAKE_CURRENT_BINARY_DIR}/../../module
                                                                   -I${CMAKE_CURRENT_BINARY_DIR}/../../module
                                                                   -I${CMAKE_CURRENT_BINARY_DIR}/../common/module
                                                                   -L${CMAKE_CURRENT_BINARY_DIR}/../../../lib
                                                                   -lcloudsc-common-lib  -l${DWARF_CLOUDSC_LIB} ${F2PY_COMMON_SOURCES} ${F2PY_DWARF_SOURCES} 
								   #> f2py_log.txt 2> f2py_log_err.txt
         COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}  ../../../bin/pythonexec/
         COMMAND ${CMAKE_COMMAND} -E copy_directory ${CMAKE_CURRENT_BINARY_DIR}/../../../lib/  ${CMAKE_CURRENT_BINARY_DIR}/../../../bin/pythonexec
         COMMAND sh -c "rm -f ../../bin/pythonexec/*.f90"
         COMMENT "Running F90wrap/F2PY after Fortran build completed"
         DEPENDS dwarf-cloudsc-lib
         VERBATIM
      )
endif()
