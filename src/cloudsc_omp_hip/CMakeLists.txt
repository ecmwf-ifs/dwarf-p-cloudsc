# (C) Copyright 1988- ECMWF.
#
# This software is licensed under the terms of the Apache Licence Version 2.0
# which can be obtained at http://www.apache.org/licenses/LICENSE-2.0.
# In applying this licence, ECMWF does not waive the privileges and immunities
# granted to it by virtue of its status as an intergovernmental organisation
# nor does it submit to any jurisdiction.

# Define this dwarf variant as an ECBuild feature
ecbuild_add_option( FEATURE CLOUDSC_OMP_HIP
	DESCRIPTION "Build the HIP version CLOUDSC using Serialbox" DEFAULT OFF
        CONDITION Serialbox_FOUND
)

if( HAVE_CLOUDSC_OMP_HIP )


    if(NOT DEFINED ROCM_PATH)
      if(DEFINED ENV{ROCM_PATH})
        set(ROCM_PATH $ENV{ROCM_PATH} CACHE PATH "Path to which ROCM has been installed")
      else()
        set(ROCM_PATH "/opt/rocm" CACHE PATH "Path to which ROCM has been installed")
      endif()
    endif()

    enable_language(HIP)
    find_package(hip REQUIRED)

    # set(CMAKE_C_COMPILER "${ROCM_PATH}/bin/hipcc")
    # set(CMAKE_CXX_COMPILER "${ROCM_PATH}/bin/hipcc")
    set(CMAKE_HIP_FLAGS "${CMAKE_HIP_FLAGS} -v -O3 -ffast-math")

    ###### SCC-HIP ####

    # cloudsc_c_k_caching.cpp  cloudsc_c_k_caching_mod.F90  CMakeLists.txt	     load_state.backup	load_state.h
    # cloudsc_c_k_caching.h	 cloudsc_driver_mod.F90       dwarf_cloudsc_gpu.F90  load_state.cpp	yoecldp_c.h

    ecbuild_add_library(
	    TARGET dwarf-cloudsc-hip-lib
	    INSTALL_HEADERS LISTED
	    SOURCES
	        load_state.h
		load_state.cpp
		yoecldp_c.h
		# cloudsc_c_k_caching.h
		cloudsc_c_k_caching.cpp
		cloudsc_c_k_caching_mod.F90
		# cloudsc_driver_mod.F90
	    PUBLIC_INCLUDES
	    	$<INSTALL_INTERFACE:include>
            	$<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/>
	    PUBLIC_LIBS
	        cloudsc-common-lib
                hip::device
		Serialbox::Serialbox_C
		$<${HAVE_OMP}:OpenMP::OpenMP_C>
		OpenMP::OpenMP_Fortran
    )

    set_source_files_properties(cloudsc_c_k_caching.cpp PROPERTIES LANGUAGE HIP)
    set_target_properties(dwarf-cloudsc-hip-lib PROPERTIES LINKER_LANGUAGE Fortran)

    # ecbuild_add_library(
    #     TARGET dwarf-cloudsc-hip-lib
    # 	INSTALL_HEADERS LISTED
    #     SOURCES
    # 	    cloudsc/yoecldp_c.h
    # 	    cloudsc/load_state.h
    #         cloudsc/load_state.cpp
    #         cloudsc/cloudsc_c.h
    #         cloudsc/cloudsc_c.cpp
    #         cloudsc/cloudsc_driver.h
    #         cloudsc/cloudsc_driver.cpp
    #         cloudsc/cloudsc_validate.h
    #         cloudsc/cloudsc_validate.cpp
    #         cloudsc/mycpu.h
    #         cloudsc/mycpu.cpp
    #     PUBLIC_INCLUDES
    #         $<INSTALL_INTERFACE:include>
    #         $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}/cloudsc>
    #     PUBLIC_LIBS
    # 	    hip::device
    #         Serialbox::Serialbox_C
    #         $<${HAVE_OMP}:OpenMP::OpenMP_C>
    # )

    target_include_directories(dwarf-cloudsc-hip-lib PUBLIC $<INSTALL_INTERFACE:include> $<BUILD_INTERFACE:${CMAKE_CURRENT_SOURCE_DIR}>)
    # target_link_libraries(dwarf-cloudsc-hip-lib PUBLIC hip::device Serialbox::Serialbox_C $<${HAVE_OMP}:OpenMP::OpenMP_C>)

    target_compile_options(dwarf-cloudsc-hip-lib PRIVATE --offload-arch=gfx90a)
    
    ecbuild_add_executable(
        TARGET  dwarf-cloudsc-omp-hip
	SOURCES dwarf_cloudsc_gpu.F90 cloudsc_driver_mod.F90
        LIBS dwarf-cloudsc-hip-lib
    )

    set_target_properties(dwarf-cloudsc-omp-hip PROPERTIES LINKER_LANGUAGE Fortran)

   # add_library(phys-hip-lib phys_hip_mod.F90 phys_hip.cpp)
   # set_source_files_properties(phys_hip.cpp PROPERTIES LANGUAGE HIP)
   # target_link_libraries(phys-hip-lib PRIVATE hip::device OpenMP::OpenMP_Fortran)
   # target_compile_options(phys-hip-lib PRIVATE "-fopenmp")
   # target_compile_options(phys-hip-lib PRIVATE "--offload-arch=gfx90a") #  -x hip")
   # add_executable( phys_LITE_LOOP_hip phys_mod.F90 phys_driver.F90 )
   # target_link_libraries(phys_LITE_LOOP_hip PRIVATE phys-hip-lib OpenMP::OpenMP_Fortran hip::device hip::host)
   # target_compile_options(phys_LITE_LOOP_hip PRIVATE "-homp") # -x hip
   # target_compile_definitions(phys_LITE_LOOP_hip PRIVATE HIP_OMP OpenMP::OpenMP_Fortran)
   # # target_link_options(phys_LITE_LOOP_hip PRIVATE "-x hip")
   # set_target_properties(phys_LITE_LOOP_hip PROPERTIES LINKER_LANGUAGE Fortran)

    
    # Create symlink for the input data
    execute_process(COMMAND ${CMAKE_COMMAND} -E create_symlink
        ${CMAKE_CURRENT_SOURCE_DIR}/../../data ${CMAKE_CURRENT_BINARY_DIR}/../../../data )

    ecbuild_add_test(
        TARGET dwarf-cloudsc-hip-serial
        COMMAND bin/dwarf-cloudsc-hip
        ARGS 1 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 1
    )
    ecbuild_add_test(
        TARGET dwarf-cloudsc-hip-omp
        COMMAND bin/dwarf-cloudsc-hip
        ARGS 4 100 16
        WORKING_DIRECTORY ${CMAKE_CURRENT_BINARY_DIR}/../../..
        OMP 4
        CONDITION HAVE_OMP
    )

endif()
