#--------------------------------------------*-cmake-*---------------------------------------------#
# file   src/CMakeLists.txt
# author Mathew Cleveland
# date   October 11th 2021
# brief  Instructions for building src level Makefile.
# note   Copyright (C) 2021 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #
cmake_minimum_required(VERSION 3.9.0)
project(odd_src_dir)

# Additional platform checks needed by Odd
include(odd-platform-checks)
odd_platform_checks()


if("${CMAKE_CXX_COMPILER_ID}" STREQUAL "Clang")
  # Suppress warnings (for now)
  list(
    APPEND
    disabled-llvm-checks
    "-Wno-enum-float-conversion"
    "-Wno-shadow"
    "-Wno-shorten-64-to-32"
    "-Wno-shadow-field"
    "-Wno-covered-switch-default"
    "-Wno-implicit-int-conversion"
    "-Wno-implicit-int-float-conversion"
    "-Wno-unused-template"
    "-Wno-header-hygiene"
    "-Wno-cast-qual"
    "-Wno-inconsistent-missing-destructor-override"
    "-Wno-missing-variable-declarations"
    "-Wno-deprecated-copy-dtor")
  foreach(flag ${disabled-llvm-checks})
    string(APPEND CMAKE_C_FLAGS " ${flag}")
    string(APPEND CMAKE_CXX_FLAGS " ${flag}")
  endforeach()
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  string(REPLACE "-Wconversion" " " CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  string(REPLACE "-Wconversion" " " CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

#
# Build components
#
message("
Adding Odd components:")

# Level 1
message(" ")
message(STATUS "Configuring Level 1 packages --")
add_dir_if_exists(odd_release)

# Level 2
message(" ")
message(STATUS "Configuring Level 2 packages --")
add_dir_if_exists(solver)

# Only process above directories in the clang-tidy CI

# Level 3
message(" ")
message(STATUS "Configuring Level 3 packages --")
add_dir_if_exists(api)

# Level 4
message(" ")
message(STATUS "Configuring Level 4 packages --")
add_dir_if_exists(driver)

# ---------------------------------------------------------------------------##
# Summary
# ---------------------------------------------------------------------------##

add_feature_info(BUILD_AUTODOC BUILD_AUTODOC
                 "Turn on to enable generation of doxygen documentation (make autodoc).")

if(VERBOSE OR ENV{VERBOSE})
  message(" ")
  feature_summary(WHAT ALL INCLUDE_QUIET_PACKAGES FATAL_ON_MISSING_REQUIRED_PACKAGES QUIET_ON_EMPTY)
endif()

if(CMAKE_CONFIGURATION_TYPES)
  set(Odd_BUILD_TYPE "Multi-config")
else()
  string(TOUPPER ${CMAKE_BUILD_TYPE} Odd_BUILD_TYPE)
endif()

message(
  "
Odd build summary:

Version     : Odd ${Odd_VERSION_FULL}
Build type  : ${Odd_BUILD_TYPE}
C4 Model    : ${DRACO_C4}
Prefix dir  : ${CMAKE_INSTALL_PREFIX}
Source dir  : ${Odd_SOURCE_DIR}
Build dir   : ${Odd_BINARY_DIR}

CXX Compiler: ${CMAKE_CXX_COMPILER}")
if(CAFS_Fortran_COMPILER)
  message("CAFS Fortran: ${CAFS_Fortran_COMPILER}")
endif()
if(CMAKE_CONFIGURATION_TYPES)
  message("CXX FLAGS   : ${CMAKE_CXX_FLAGS} ")
  message("CXX Debug FL: ${CMAKE_CXX_FLAGS_DEBUG}")
  message("CXX Release : ${CMAKE_CXX_FLAGS_RELEASE}")
  if(_LANGUAGES_ MATCHES Fortran)
    message("Fortran     : ${CMAKE_Fortran_COMPILER}")
    message("FC FLAGS    : ${CMAKE_Fortran_FLAGS}")
    message("FC Debug FLA: ${CMAKE_Fortran_FLAGS_DEBUG}")
    message("FC Release F: ${CMAKE_Fortran_FLAGS_RELEASE}")
  endif()
  if(_LANGUAGES_ MATCHES CUDA)
    message("CUDA        : ${CMAKE_CUDA_COMPILER}")
    message("Cuda FLAGS  : ${CMAKE_CUDA_FLAGS}")
    message("Cuda Debug F: ${CMAKE_CUDA_FLAGS_DEBUG}")
    message("Cuda Release: ${CMAKE_CUDA_FLAGS_RELEASE}")
  endif()
else()
  message("C FLAGS     : ${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${Odd_BUILD_TYPE}}")
  message("CXX FLAGS   : ${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${Odd_BUILD_TYPE}}")
  if(_LANGUAGES_ MATCHES Fortran)
    message("Fortran     : ${CMAKE_Fortran_COMPILER}")
    message("FC FLAGS    : ${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${Odd_BUILD_TYPE}}")
  endif()
  if(_LANGUAGES_ MATCHES CUDA)
    message("CUDA        : ${CMAKE_CUDA_COMPILER} (${CMAKE_CUDA_COMPILER_ID})")
    message("Cuda FLAGS  : ${CMAKE_CUDA_FLAGS} ${CMAKE_CUDA_FLAGS_${Draco_BUILD_TYPE}}")
  endif()
endif()
if(CAFS_Fortran_COMPILER)
  message("CAFS Fortran: ${CAFS_Fortran_COMPILER}")
endif()
if("${DRACO_C4}" STREQUAL "MPI")
  message("mpirun cmd  : ${MPIEXEC_EXECUTABLE} ${MPIEXEC_NUMPROC_FLAG} N ${MPIEXEC_PREFLAGS}")
endif()
message(
  "
Library Type: ${DRACO_LIBRARY_TYPE}

Draco       : ${DRACO_DIR}
")
if(CMAKE_CXX_COMPILER_WRAPPER STREQUAL CrayPrgEnv)
  message("Cray system detected: CC -craype-verbose -V|--version:
")
  if(${CMAKE_CXX_COMPILER_ID} STREQUAL "Cray" AND NOT DEFINED ENV{CRAY_PE_USE_CLANG})
    set(ver_opt "-V")
  else()
    set(ver_opt "--version")
  endif()
  execute_process(COMMAND ${CMAKE_CXX_COMPILER} -craype-verbose ${ver_opt} OUTPUT_VARIABLE tmp)
  message("${tmp}")
endif()

# Push some variables up one level
set(Odd_EXPORT_TARGET_PROPERTIES
    "${Odd_EXPORT_TARGET_PROPERTIES}"
    PARENT_SCOPE)

# ------------------------------------------------------------------------------------------------ #
# End of src/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
