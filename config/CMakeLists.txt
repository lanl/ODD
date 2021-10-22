#--------------------------------------------*-cmake-*---------------------------------------------#
# file   config/CMakeLists.txt
# author Mathew Cleveland
# date   October 21st 2021
# brief  Install the config directory in the target directory
# note   Copyright (C) 2010-2021 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #
cmake_minimum_required(VERSION 3.14.0)
project(config)

include(CMakePackageConfigHelpers)

file(GLOB CMake_src *.cmake)
file(GLOB CMake_in *.in)

string(TOUPPER "${CMAKE_BUILD_TYPE}" upper_build_type)
set(Odd_CXX_COMPILER_FLAGS "${CMAKE_CXX_FLAGS} ${CMAKE_CXX_FLAGS_${upper_build_type}}")
set(Odd_C_COMPILER_FLAGS "${CMAKE_C_FLAGS} ${CMAKE_C_FLAGS_${upper_build_type}}")
set(Odd_Fortran_COMPILER_FLAGS
    "${CMAKE_Fortran_FLAGS} ${CMAKE_Fortran_FLAGS_${upper_build_type}}")

# Publish the active build_types (A subset of CMAKE_CONFIGURATION_TYPES)
if(CMAKE_CONFIGURATION_TYPES)
 unset(ODD_CONFIGURATION_TYPES CACHE)
  foreach(cfg ${CMAKE_CONFIGURATION_TYPES})
    if(EXISTS ${Odd_BINARY_DIR}/${cfg})
      if("${ODD_CONFIGURATION_TYPES}x" STREQUAL "x")
        set(ODD_CONFIGURATION_TYPES "${cfg}")
      else()
        list(APPEND ODD_CONFIGURATION_TYPES ${cfg})
      endif()
    endif()
  endforeach()
  set(ODD_CONFIGURATION_TYPES
      "${ODD_CONFIGURATION_TYPES}"
      CACHE STRING "Available multiconfig builds.")
endif()

# Process the config file
configure_file(odd-config-install.cmake.in "${CMAKE_CURRENT_BINARY_DIR}//odd-config.cmake"
               @ONLY)
write_basic_package_version_file(
  "${CMAKE_CURRENT_BINARY_DIR}/odd-config-version.cmake"
  VERSION ${Odd_VERSION}
  COMPATIBILITY AnyNewerVersion)

# Install scripts and macros to make them available by other projects.
if(CMAKE_CONFIGURATION_TYPES)
  set(file_list ${CMake_src} ${CMake_in})
else()
  set(file_list ${CMake_src} ${CMake_in} "${CMAKE_CURRENT_BINARY_DIR}/odd-config.cmake"
                "${CMAKE_CURRENT_BINARY_DIR}/odd-config-version.cmake")
endif()

install(FILES ${file_list} DESTINATION cmake)

# ------------------------------------------------------------------------------------------------ #
# End config/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
