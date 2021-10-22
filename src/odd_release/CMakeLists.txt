#--------------------------------------------*-cmake-*---------------------------------------------#
# file   odd_release/CMakeLists.txt
# author Kelly Thompson <kgt@lanl.gov>
# date   2012 Aug 1
# brief  Generate build project files for odd_release.
# note   Copyright (C) 2015-2021 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #
cmake_minimum_required(VERSION 3.9.0)
project(odd_release CXX)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  string(REPLACE "/W2" "/W4" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  string(REPLACE "/W2" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  string(APPEND CMAKE_C_FLAGS " -Wconversion")
  string(APPEND CMAKE_CXX_FLAGS " -Wconversion")
endif()

# ------------------------------------------------------------------------------------------------ #
# Generate config.h (only occurs when cmake is run)
# ------------------------------------------------------------------------------------------------ #

string(TIMESTAMP Odd_BUILD_DATE "%Y %b %d")
configure_file(config.h.in ${PROJECT_BINARY_DIR}/odd_release/config.h)

# Generate the dll_declspec.h file that defines DLL_PUBLIC_<pkg> CPP macros.
include(generate_dll_declspec)
set(dirs_with_libs
    odd_release
    solver
    driver)
generate_dll_declspec("odd_release" "${dirs_with_libs}")

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #
file(GLOB sources Release.cc)
file(GLOB headers Release.hh)
list(APPEND headers ${PROJECT_BINARY_DIR}/odd_release/config.h
     ${PROJECT_BINARY_DIR}/odd_release/dll_declspec.h)

# ------------------------------------------------------------------------------------------------ #
# Build package library
# ------------------------------------------------------------------------------------------------ #
set(include_dirs
    PUBLIC $<BUILD_INTERFACE:${Odd_SOURCE_DIR}/src>
    $<BUILD_INTERFACE:${PROJECT_BINARY_DIR}> # config.h
    $<INSTALL_INTERFACE:${DBSCFG_IMPORT_PREFIX}>)

add_component_library(
  TARGET Lib_odd_release
  TARGET_DEPS "Lib_diagnostics"
  INCLUDE_DIRS "${include_dirs}"
  LIBRARY_NAME odd_release
  SOURCES "${sources}"
  HEADERS "${headers}"
  EXPORT_NAME "odd-targets")

# Copy necessary Draco dll files to the CMAKE_RUNTIME_OUTPUT_DIRECTORY directory
if(MSVC AND DRACO_LIBRARY_TYPE STREQUAL SHARED)
  foreach(draco_dll Lib_c4 Lib_dsxx Lib_diagnostics)
    add_custom_command(
      TARGET Lib_odd_release
      POST_BUILD
      COMMAND ${CMAKE_COMMAND} -E copy_if_different $<TARGET_FILE:${draco_dll}>
              $<TARGET_FILE_DIR:Lib_odd_release>
      COMMENT "Copying dll files to output directory to allow unit tests to run.")
  endforeach()
endif()

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #
install(FILES ${headers} DESTINATION include/odd_release)

# ------------------------------------------------------------------------------------------------ #
# Unit tests
# ------------------------------------------------------------------------------------------------ #
if(BUILD_TESTING)
  add_subdirectory(test)
endif()

# ------------------------------------------------------------------------------------------------ #
# End odd_release/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
