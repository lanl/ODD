# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   odd/CMakeLists.txt
# author Mathew Cleveland (Pulled from the Jayenne project)
# date   October 11th, 2021
# brief  Instructions for building root level Makefile.
# note   Copyright (C) 2013-2021 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #
cmake_minimum_required(VERSION 3.17.0 FATAL_ERROR)

string(CONCAT odd_descript "This project provides a simple TRT orthogonal diffusion discretization (odd) solver interface library. ")
project(
  Odd
  VERSION 0.0
  DESCRIPTION ${odd_descript}
  LANGUAGES CXX C)

# Don't warn about "An imported target missing its location property". Only for Cmake-3.19+.
# cmake_policy(SET CMP0111 OLD)

# Create '${PROJECT_BINARY_DIR}/compile_commands.json'.  Used by clang-analyze and clang-tidy.
set(CMAKE_EXPORT_COMPILE_COMMANDS TRUE)

# Build system configuration files are located here
set(CMAKE_MODULE_PATH ${PROJECT_SOURCE_DIR}/config)

# Provide default value for install_prefix
if("${CMAKE_INSTALL_PREFIX}" STREQUAL "/usr/local" OR "${CMAKE_INSTALL_PREFIX}" MATCHES
                                                      "C:/Program Files")
  set(CMAKE_INSTALL_PREFIX "${CMAKE_CURRENT_BINARY_DIR}/../install")
  get_filename_component(CMAKE_INSTALL_PREFIX "${CMAKE_INSTALL_PREFIX}" ABSOLUTE)
  set(CMAKE_INSTALL_PREFIX
      "${CMAKE_INSTALL_PREFIX}"
      CACHE PATH "Install path prefix, prepended onto install directories" FORCE)
endif()

#
# DRACO
#
# Locate a Draco installation (find and source 'draco-config.cmake')
foreach(dpath ${CMAKE_INSTALL_PREFIX} $ENV{DRACO_DIR} ${DRACO_DIR})
  file(TO_CMAKE_PATH "${dpath}" draco_config_dir)
  if(EXISTS ${draco_config_dir})
    # build system files are located here.  This allows include(dracoVersion), etc.
    list(APPEND CMAKE_MODULE_PATH ${draco_config_dir})
    # find_package(draco) needs this
    list(APPEND CMAKE_PREFIX_PATH ${draco_config_dir})
  endif()
  unset(draco_config_dir)
endforeach()
# This command sets the variable draco_DIR
find_package(draco 7.11 REQUIRED CONFIG)
get_filename_component(DRACO_DIR ${draco_DIR} DIRECTORY)

# Now that Draco's Lib_dsxx is known to Odd, build a list of regular expressions that gcov
# should ignore. (These commands must be located in CMakeLists.txt after draco is discovered, but
# before 'dbs_basic_setup' is called.  The make target 'covrep' is created when dbs_basic_setup ->
# compilerEnv.cmake and will use the current value of CODE_COVERAGE_IGNORE_REGEX).
if(CODE_COVERAGE
   AND UNIX
   AND CMAKE_CXX_COMPILER_ID STREQUAL GNU)
  # Exclude installed Draco headers from coverage report.
  get_target_property(cov_excl_dir Lib_dsxx INTERFACE_INCLUDE_DIRECTORIES)
  if(EXISTS ${cov_excl_dir})
    list(
      APPEND
      CODE_COVERAGE_IGNORE_REGEX
      ${cov_excl_dir}/*
      *config.h)
  endif()
  unset(cov_excl_dir)
endif()

# Setup build environment by calling macros from draco-config.cmake. The following Draco cmake files
# are sourced: - platform_checks.cmake - compilerEnv.cmake - sets compiler flags (C, CXX, Fortran,
# Cuda) and creates 'covrep' target. - component_macros.cmake - vendor_libraries (find third party
# libraries).
dbs_basic_setup()

#
# Odd version
#
include(dracoVersion)
set_ccs2_software_version(Odd)

#
# Unit Test Setup
#
include(dracoTesting)
if(RUN_VERIFICATION_TESTS)
  set(CTEST_TESTING_TIMEOUT "7200") # 2 hours
else()
  set(CTEST_TESTING_TIMEOUT "300") # 5 minutes
endif()

# Set default compile environment:

# Setup defaults, value checks, etc.
include(buildEnv)
# set defaults for BUILD_TYPE and INSTALL_PREFIX
dbssetdefaults()
# Initialize data for exporting Odd targets
dbsinitexporttargets("Odd")
# Save config info
dbsconfiginfo()

# Extra compiler setup
dbssetupstaticanalyzers()

#
# Build the components:
#
add_subdirectory(src)
add_subdirectory(autodoc)
add_subdirectory(config)

# install top level documents
install(FILES ChangeLog LICENSE.md README.md DESTINATION ${CMAKE_INSTALL_PREFIX})

# Export targets
install(EXPORT odd-targets DESTINATION cmake)

# If the build system is building all runtime objects (dll and exe files) in a common directory, add
# a rule to copy all of the Draco dll files to CMAKE_RUNTIME_OUTPUT_DIRECTORY
if(DEFINED CMAKE_RUNTIME_OUTPUT_DIRECTORY)
  unset(dll_copy_command)
  foreach(dlib ${Draco_LIBRARIES})
    list(APPEND dll_copy_command COMMAND
         ${CMAKE_COMMAND};-E;copy_if_different;$<TARGET_FILE:${dlib}>;$<TARGET_FILE_DIR:Lib_mc>)
  endforeach()
  add_custom_target(
    copy_draco_dll ALL
    ${dll_copy_command}
    COMMENT "Copying dll files to CMAKE_RUNTIME_OUTPUT_DIRECTORY.")
endif()

#
# Extra custom commands
#

# Build target: "vtest" Reconfigure the build directory to include the long running tests and run
# the tests.
if(RUN_VERIFICATION_TESTS)
  add_custom_target(
    vtest
    COMMAND "${CMAKE_CTEST_COMMAND}" -L vtest --force-new-ctest-process ${ARGS}
    COMMENT "Running verification tests (label 'vtest')")
else()
  add_custom_target(
    vtest
    # reconfigure with RUN_VERIFICATION_TESTS=ON
    COMMAND "${CMAKE_COMMAND}" -DRUN_VERIFICATION_TESTS=ON "${PROJECT_SOURCE_DIR}"
    # rebuild with RUN_VERIFICATION_TESTS=ON
    COMMAND "${CMAKE_COMMAND}" --build "${PROJECT_BINARY_DIR}"
    # run the perfbench tests
    COMMAND "${CMAKE_CTEST_COMMAND}" -L vtest --force-new-ctest-process ${ARGS}
    COMMENT "Enabling long running verification test mode and reconfiguring...")
endif()
add_feature_info(RUN_VERIFICATION_TESTS RUN_VERIFICATION_TESTS
                 "Run long running verification tests along with normal unit tests.")
# ------------------------------------------------------------------------------------------------ #
# End CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
