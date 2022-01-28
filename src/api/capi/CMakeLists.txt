# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   api/CMakeLists.txt
# author Mathew Cleveland
# date   October 2021
# brief  Generate build project files for mc
# note   Copyright (C) 2015-2021 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #
cmake_minimum_required(VERSION 3.9.0)
project(capi C)

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

file(GLOB headers *.h)
file(GLOB sources *.c)

if(CMAKE_CXX_COMPILER_ID STREQUAL "MSVC")
  string(REPLACE "/W2" "/W4" CMAKE_C_FLAGS "${CMAKE_C_FLAGS}")
  string(REPLACE "/W2" "/W4" CMAKE_CXX_FLAGS "${CMAKE_CXX_FLAGS}")
endif()

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #
install(FILES ${headers} DESTINATION include/capi)

# ------------------------------------------------------------------------------------------------ #
# Unit tests
# ------------------------------------------------------------------------------------------------ #
if(BUILD_TESTING)
 add_subdirectory(test)
endif()

# ------------------------------------------------------------------------------------------------ #
# Autodoc
# ------------------------------------------------------------------------------------------------ #
process_autodoc_pages()

# ------------------------------------------------------------------------------------------------ #
# End mc/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
