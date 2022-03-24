#--------------------------------------------*-cmake-*---------------------------------------------#
# file   solver/test/CMakeLists.txt
# author Mathew Cleveland
# date   October 21st 2021
# brief  Generate build project files for solver/test.
# note   Copyright (C) 2015-2021 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #
project(solver_test CXX)

# ------------------------------------------------------------------------------------------------ #
# Source files
# ------------------------------------------------------------------------------------------------ #

file(GLOB test_sources *.cc)

# ------------------------------------------------------------------------------------------------ #
# Build Unit tests
# ------------------------------------------------------------------------------------------------ #


set(two_proc_tests
    ${PROJECT_SOURCE_DIR}/tstGhost_Comm.cc)
list(REMOVE_ITEM test_sources ${two_proc_tests})

# Remaining tests will be run with 1, 2 and 4 PEs
add_parallel_tests(
  SOURCES "${test_sources}"
  PE_LIST "1;2;4"
  DEPS "Lib_solver")

add_parallel_tests(
  SOURCES "${two_proc_tests}"
  PE_LIST "2"
  DEPS "Lib_solver")


file(GLOB ipcress_files ${PROJECT_SOURCE_DIR}/*.ipcress)
provide_aux_files(FILES "${ipcress_files}" FOLDER solver_test)


# ------------------------------------------------------------------------------------------------ #
# end CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
