# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   driver/CMakeLists.txt
# author Mathew Cleveland
# date   October 2021
# brief  Generate build project files for mc
# note   Copyright (C) 2015-2021 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #
cmake_minimum_required(VERSION 3.9.0)
project(driver CXX)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  string(REPLACE "/W2" "/W4" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  string(REPLACE "/W2" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
elseif("${CMAKE_CXX_COMPILER_ID}" STREQUAL "GNU")
  string(APPEND CMAKE_C_FLAGS " -Wconversion")
  string(APPEND CMAKE_CXX_FLAGS " -Wconversion")
endif()

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #

file(GLOB headers *.hh)
file(GLOB template_implementations *.t.hh *.i.hh)
file(GLOB sources *.cc)
file(GLOB explicit_instantiations *_pt.cc)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  string(REPLACE "/W2" "/W4" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  string(REPLACE "/W2" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

# ------------------------------------------------------------------------------------------------ #
# Build package library
# ------------------------------------------------------------------------------------------------ #
add_component_library(
  TARGET Lib_driver
  TARGET_DEPS "Lib_solver"
  LIBRARY_NAME driver
  SOURCES "${sources}"
  HEADERS "${headers}"
  EXPORT_NAME "odd-targets")

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #
install(FILES ${headers} DESTINATION include/driver)

# ------------------------------------------------------------------------------------------------ #
# Unit tests
# ------------------------------------------------------------------------------------------------ #
## if(BUILD_TESTING)
##   add_subdirectory(test)
## endif()
## 
# ------------------------------------------------------------------------------------------------ #
# Autodoc
# ------------------------------------------------------------------------------------------------ #
process_autodoc_pages()

# ------------------------------------------------------------------------------------------------ #
# End mc/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
