# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

### Helper macro to create a target and associated tests

macro( cloudsc_add_cuda_exe TARGET DRIVER KERNEL )

    ecbuild_add_library(
        TARGET ${TARGET}-${prec}-lib
        INSTALL_HEADERS LISTED
        SOURCES
            cloudsc/yoecldp_c.h
            cloudsc/dtype.h
            cloudsc/load_state.h
            cloudsc/load_state.cu
            cloudsc/cloudsc_validate.h
            cloudsc/cloudsc_validate.cu
            cloudsc/mycpu.h
            cloudsc/mycpu.cu
            ${DRIVER}
            ${KERNEL}
        PUBLIC_INCLUDES
            $<INSTALL_INTERFACE:include>
            $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>
        PUBLIC_LIBS
            $<${HAVE_HDF5}:hdf5::hdf5>
            $<${HAVE_SERIALBOX}:Serialbox::Serialbox_C>
            $<${HAVE_OMP}:OpenMP::OpenMP_C>
            CUDA::cudart
        DEFINITIONS
            ${CLOUDSC_DEFINITIONS}
            $<$<STREQUAL:"${prec}","sp">:SINGLE>
    )

    target_compile_options( ${TARGET}-${prec}-lib
        PRIVATE $<$<COMPILE_LANGUAGE:CUDA>:SHELL:${CLOUDSC_CUDA_OPT_FLAGS} ${CLOUDSC_CUDA_FLAGS}>
    )
    set_target_properties( ${TARGET}-${prec}-lib
        PROPERTIES CUDA_SEPARABLE_COMPILATION ON
    )

    ecbuild_add_executable(
        TARGET      ${TARGET}-${prec}
        SOURCES     dwarf_cloudsc.cpp
        LIBS        ${TARGET}-${prec}-lib
    )

    ecbuild_add_test(
        TARGET      ${TARGET}-${prec}-serial
        COMMAND     ${TARGET}-${prec}
        ARGS        1 1000 128
        WORKING_DIRECTORY ${CMAKE_BINARY_DIR}
    )

endmacro()

### Kernel variants

# SCC
cloudsc_add_cuda_exe(
    dwarf-cloudsc-gpu-cuda-c
    "cloudsc/cloudsc_driver.h;cloudsc/cloudsc_driver.cu"
    "cloudsc/cloudsc_c.h;cloudsc/cloudsc_c.cu"
)

# SCC-HOIST
cloudsc_add_cuda_exe(
    dwarf-cloudsc-gpu-cuda-c-hoist
    "cloudsc/cloudsc_driver_hoist.h;cloudsc/cloudsc_driver_hoist.cu"
    "cloudsc/cloudsc_c_hoist.h;cloudsc/cloudsc_c_hoist.cu"
)

# SCC-K-CACHING
cloudsc_add_cuda_exe(
    dwarf-cloudsc-gpu-cuda-c-k-caching
    "cloudsc/cloudsc_driver.h;cloudsc/cloudsc_driver.cu"
    "cloudsc/cloudsc_c_k_caching.h;cloudsc/cloudsc_c_k_caching.cu"
)

# SCC-OPT
cloudsc_add_cuda_exe(
    dwarf-cloudsc-gpu-cuda-c-opt
    "cloudsc/cloudsc_driver_opt.h;cloudsc/cloudsc_driver_opt.cu"
    "cloudsc/cloudsc_c_opt.h;cloudsc/cloudsc_c_opt.cu"
)
