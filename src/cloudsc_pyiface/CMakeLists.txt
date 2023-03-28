# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_FORTRAN_PYIFACE
    DESCRIPTION "Build the Python interface to the (slightly customized) cloudsc-fortran" DEFAULT OFF
    CONDITION HDF5_FOUND
)
ecbuild_add_option( FEATURE CLOUDSC_FORTRAN_PYIFACE_BINARY
    DESCRIPTION "Build the debug-oriented binary for Python interface of (customized) cloudsc-fortran" DEFAULT OFF
    CONDITION HAVE_CLOUDSC_FORTRAN_PYIFACE
)

if( HAVE_CLOUDSC_FORTRAN_PYIFACE )

    # Utilities to manage Python virtual environments
    include( python_venv )

    # Set up a custom venv for this variant and install the necessary dependencies
    set( pyiface_VENV_PATH ${CMAKE_BINARY_DIR}/venv_pyiface )
    setup_python_venv( ${pyiface_VENV_PATH} )

    if( NOT Python3_EXECUTABLE )
        ecbuild_error("[PyIface] Could not find Python3 executable in virtualenv")
    endif()

    # Install the f90wrap build dependency and via latest pip
    execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install --upgrade pip)
    execute_process(COMMAND ${Python3_EXECUTABLE} -m pip install -e ${CMAKE_CURRENT_SOURCE_DIR})

    # Define module directory to facilitate f90wrap/f2py execution
    ecbuild_enable_fortran(MODULE_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../module)
    message (STATUS "Module directory set to:")
    message (STATUS  ${CMAKE_Fortran_MODULE_DIRECTORY})


    # Set names and lists to abstract the f90wrap/f2py process
    set( DWARF_CLOUDSC_LIB  dwarf-cloudsc-lib)
    set( DWARF_COMMON_LIB  cloudsc-common-lib)
    set( PYTHON_MODN cloudsc)
    set( FORTRAN_PYTHON_COMMON_SOURCES
        yomphyder yoecldp yoethf yomcst yoephli expand_mod
    )
    set( FORTRAN_PYTHON_DWARF_SOURCES cloudsc_driver_mod )
    set( COMMON_MOD_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/../common/module/ )
    set( COMMON_MOD_BIN_LOCATION ${CMAKE_CURRENT_BINARY_DIR}/../common/module/ )
    set( FORTRAN_SRC_LOCATION ${CMAKE_CURRENT_SOURCE_DIR}/fortransrc/ )

    # Manipulate lists to prepare abstracted command-line input/output for f90wrap/f2py
    set( F90WRAP_COMMON_SOURCES ${FORTRAN_PYTHON_COMMON_SOURCES})
    set(  F90WRAP_DWARF_SOURCES  ${FORTRAN_PYTHON_DWARF_SOURCES})
    set(    F2PY_COMMON_SOURCES ${FORTRAN_PYTHON_COMMON_SOURCES})
    set(     F2PY_DWARF_SOURCES  ${FORTRAN_PYTHON_DWARF_SOURCES})
    list(TRANSFORM F90WRAP_COMMON_SOURCES PREPEND ${COMMON_MOD_LOCATION})
    list(TRANSFORM F90WRAP_COMMON_SOURCES APPEND .F90 )
    list(TRANSFORM F90WRAP_DWARF_SOURCES  PREPEND ${FORTRAN_SRC_LOCATION})
    list(TRANSFORM F90WRAP_DWARF_SOURCES  APPEND .F90 )
    list(TRANSFORM F2PY_COMMON_SOURCES    PREPEND f90wrap_)
    list(TRANSFORM F2PY_COMMON_SOURCES    APPEND .f90 )
    list(TRANSFORM F2PY_DWARF_SOURCES     PREPEND f90wrap_)
    list(TRANSFORM F2PY_DWARF_SOURCES     APPEND .f90 )

    # Build CLOUDSC driver/kernel library, to be further linked by f2py
    ecbuild_add_library( TARGET ${DWARF_CLOUDSC_LIB}
        SOURCES
            ./fortransrc/cloudsc_driver_mod.F90
            ./fortransrc/cloudsc.F90
        PUBLIC_LIBS
             ${DWARF_COMMON_LIB}
        DEFINITIONS ${CLOUDSC_DEFINITIONS}
    )

    # Two-step F90wrap process: First, generate Fortran wrappers
    add_custom_command( TARGET ${DWARF_CLOUDSC_LIB} POST_BUILD
        COMMAND ${Python3_VENV_BIN}/f90wrap -m${PYTHON_MODN}
                    ${F90WRAP_COMMON_SOURCES} ${F90WRAP_DWARF_SOURCES}
                    -k ${CMAKE_CURRENT_SOURCE_DIR}/kind_map
                    > f90wrap_log.txt 2> f90wrap_log_err.txt
        COMMENT "[PyIface] Executing f90wrap to generate Fortran wrappers"
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        VERBATIM
    )

    # Two-step F90wrap process: Then compile and generate Python wrappers via F2Py
    # Note that we execute this in CMAKE_BINARY_DIR to make the resulting
    # _cloudsc.arch.so library directly available for dynamic loading, as f2py-f90wrap
    # does not allow specifying a particular output directory or path.
    add_custom_command(TARGET ${DWARF_CLOUDSC_LIB}  POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E env --unset=LD_LIBRARY_FLAGS
                    LDFLAGS=-Wl,-rpath,${CMAKE_BINARY_DIR}/lib
                    NPY_DISTUTILS_APPEND_FLAGS=1
                    ${Python3_VENV_BIN}/f2py-f90wrap -c
                    --f90exec=${CMAKE_Fortran_COMPILER}
                    --f90flags=${CMAKE_Fortran_FLAGS}
                    --f90flags=${CMAKE_Fortran_FLAGS_${CMAKE_BUILD_TYPE_CAPS}}
                    -m _${PYTHON_MODN}
                    -I${CMAKE_Fortran_MODULE_DIRECTORY}
                    -I${COMMON_MOD_BIN_LOCATION}
                    -L${CMAKE_BINARY_DIR}/lib
                    -l${DWARF_COMMON_LIB}
                    -l${DWARF_CLOUDSC_LIB}
                    ${F2PY_COMMON_SOURCES} ${F2PY_DWARF_SOURCES}
                    > f2py_log.txt 2> f2py_log_err.txt
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        COMMENT "[PyIface] Executing f2py-f90wrap to compile and generate Python wrappers"
        VERBATIM
    )

    # Copy the CLI driver script into the bin directory for execution
    add_custom_command(TARGET ${DWARF_CLOUDSC_LIB} POST_BUILD
        COMMAND ${CMAKE_COMMAND} -E make_directory ${CMAKE_BINARY_DIR}/bin
        COMMAND ${CMAKE_COMMAND} -E create_symlink ${Python3_VENV_BIN}/cloudsc_pyiface.py ${CMAKE_BINARY_DIR}/bin/cloudsc_pyiface.py
        COMMENT "[PyIface] Installing Python package and driver via 'pip install'"
    )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-pyiface
        COMMAND bin/cloudsc_pyiface.py
        ARGS --numomp=1 --ngptot=100 --nproma=16 --cloudsc-path=${CMAKE_BINARY_DIR} --input-path=${CMAKE_BINARY_DIR}
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
        OMP 1
    )


    if( HAVE_CLOUDSC_FORTRAN_PYIFACE_BINARY )
        # Define the (optional) binary build target for this variant
        ecbuild_add_executable( TARGET dwarf-cloudsc-fortran-pyiref
            SOURCES
                ./fortransrc/dwarf_cloudsc.F90
            LIBS
                ${DWARF_COMMON_LIB}
                ${DWARF_CLOUDSC_LIB}
            DEFINITIONS ${CLOUDSC_DEFINITIONS}
         )

        ecbuild_add_test(
            TARGET dwarf-cloudsc-fortran-pyiref
            COMMAND bin/dwarf-cloudsc-fortran-pyiref
            ARGS 1 100 16
            WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
            OMP 1
        )
    endif()

    if( HAVE_HDF5 )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/input.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../input.h5 )
      execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
          ${CMAKE_CURRENT_SOURCE_DIR}/../../config-files/reference.h5 ${CMAKE_CURRENT_BINARY_DIR}/../../../reference.h5 )
    endif()

endif()
