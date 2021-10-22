# ------------------------------------------------------------------------------------------------ #
# file  : add_driver_test.cmake
# author: Kelly Thompson <kgt@lanl.gov>
# date  : 2011
# brief : Register a driver-based test (mcgrid or ncpd)
# note  : Copyright (C) 2011-2021 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

# ------------------------------------------------------------------------------------------------ #
# add_driver_test
# ------------------------------------------------------------------------------------------------ #

# cmake-lint: disable=R0912,R0915
macro(add_driver_test input)

  # Keep a global count of how many times this macro is called.
  if(NOT DEFINED orte_tmpdir_base_enum)
    set(orte_tmpdir_base_enum 0)
  else()
    math(EXPR orte_tmpdir_base_enum "${orte_tmpdir_base_enum} + 1")
  endif()
  set(orte_tmpdir_base_enum
      ${orte_tmpdir_base_enum}
      CACHE INTERNAL "help openmpi")

  cmake_parse_arguments(
    amt "CHECK_RUN_TIME;VERIFY_MODE"
    "EPSILON;INPUT;BENCHFILE;EXEC;PROC;RESTART_DIR;RESOURCE_LOCK;TIMEOUT"
    "NUM_PROC;RUN_AFTER;LABEL" ${ARGV})

  # Enable runtime checks
  if(${amt_CHECK_RUN_TIME})
    set(CHECK_RUN_TIME "-t")
  endif()

  # Default verification mode to false if not found in arg list
  if(${amt_VERIFY_MODE})
    set(VERIFY_MODE "-y")
  endif()

  # Tolerance for comparing output.
  if(${amt_EPSILON})
    set(EPSILON "-e ${amt_EPSILON}")
  endif()

  # Format resource lock command
  if(NOT "${amt_RESOURCE_LOCK}none" STREQUAL "none")
    set(amt_RESOURCE_LOCK "RESOURCE_LOCK ${amt_RESOURCE_LOCK}")
  endif()

  # Due to naming conventions, we can deduce the name and location of the appropriate mcgrid binary:
  string(REPLACE "_test" "" mcgrid_dir ${PROJECT_NAME})
  # We no longer have mesh dependent versions

  # If the file specified by amt_INPUT is not available, try other paths:
  if(NOT EXISTS ${amt_INPUT} AND EXISTS ${PROJECT_SOURCE_DIR}/${amt_INPUT})
    set(amt_INPUT ${PROJECT_SOURCE_DIR}/${amt_INPUT})
  endif()
  if(NOT EXISTS ${amt_INPUT} AND EXISTS ${PROJECT_BINARY_DIR}/${amt_INPUT})
    set(amt_INPUT ${PROJECT_BINARY_DIR}/${amt_INPUT})
  endif()

  # if the file specified by amt_BENCHFILE is not available, try other paths:
  if(NOT EXISTS ${amt_BENCHFILE} AND EXISTS ${PROJECT_SOURCE_DIR}/${amt_BENCHFILE})
    set(amt_BENCHFILE ${PROJECT_SOURCE_DIR}/${amt_BENCHFILE})
  endif()

  # for historical reasons, default num_proc = 2
  if(NOT DEFINED amt_NUM_PROC)
    set(amt_NUM_PROC 2)
  endif()

  # short name for test
  get_filename_component(raw_name "${amt_INPUT}" NAME)

  # trim the py extension for the name registration
  string(REPLACE ".py" "" short_name ${raw_name})

  # special name for restart_dir
  if(NOT DEFINED amt_RESTART_DIR)
    set(amt_RESTART_DIR ${PROJECT_BINARY_DIR}/${short_name}_restart)
  endif()

  # Parallel or scalar
  if("${DRACO_C4}" STREQUAL "SCALAR")
    set(run_scalar "-s")
    set(amt_NUM_PROC 1)
  endif()

  # make sure python and run_driver_test.py are available
  if(NOT EXISTS ${Odd_BINARY_DIR}/tools/run_driver_test.py)
    message(FATAL_ERROR "Could not find the python Mcgrid test driver script:"
                        " ${Odd_BINARY_DIR}/tools/run_driver_test.py")
  endif()
  if(NOT EXISTS ${Python_EXECUTABLE})
    message(FATAL_ERROR "Could not find python.")
  endif()

  # register the test, if multiple processors were specied register additional tests
  list(LENGTH amt_NUM_PROC n_tests)

  foreach(iproc ${amt_NUM_PROC})
    set(multi_short_name ${short_name})
    if(n_tests GREATER 1)
      set(multi_short_name "${short_name}_${iproc}")
    endif()
    add_test(
      NAME ${mcgrid_dir}_${multi_short_name}
      COMMAND
        ${Python_EXECUTABLE} ${Odd_BINARY_DIR}/tools/run_driver_test.py
        # -x $<TARGET_FILE:Exe_mcgrid>
        -x ${amt_EXEC} -p ${iproc} -i ${amt_INPUT} # input file
        -w ${PROJECT_BINARY_DIR} # work dir
        -b ${amt_BENCHFILE} # gold standard file
        -r ${amt_RESTART_DIR} -q ${orte_tmpdir_base_enum} ${CHECK_RUN_TIME} # blank, or '-t' for
                                                                            # perfbench
        ${VERIFY_MODE} # blank, or '-y' for verification-mode
        ${EPSILON}
        # -u           # Updates gold standard file (careful!!) You will need to run 'make
        # rebuild_cache' after modifying this option.
        ${run_scalar} # run scalar (-s) or mpirun (default)
    )

    set(num_procs ${iproc})

    set_tests_properties(
      ${mcgrid_dir}_${multi_short_name}
      PROPERTIES PASS_REGULAR_EXPRESSION
                 ".*[Tt]est PASSED.*"
                 FAIL_REGULAR_EXPRESSION
                 ".*[Tt]est FAILED.*"
                 WORKING_DIRECTORY
                 "${PROJECT_BINARY_DIR}"
                 PROCESSORS
                 ${iproc}
                 # RUN_SERIAL              ON
                 LABELS
                 "nomemcheck;${amt_LABEL}")

    if(DEFINED amt_RUN_AFTER)
      string(REPLACE ".py" "" short_name_after ${amt_RUN_AFTER})
      set_tests_properties(${mcgrid_dir}_${multi_short_name}
                           PROPERTIES DEPENDS ${mcgrid_dir}_${short_name_after})
    endif()

    if(DEFINED amt_RESOURCE_LOCK)
      set_tests_properties(${mcgrid_dir}_${multi_short_name} PROPERTIES RESOURCE_LOCK
                                                                        "${amt_RESOURCE_LOCK}")
    endif()

    if(DEFINED amt_TIMEOUT)
      set_tests_properties(${mcgrid_dir}_${multi_short_name} PROPERTIES TIMEOUT "${amt_TIMEOUT}")
    endif()
  endforeach()

endmacro()

# Wrapper for adding a mcgrid test (python driver)
macro(add_mcgrid_test input)
  add_driver_test(${ARGV} EXEC $<TARGET_FILE:Exe_mcgrid>)
endmacro()

# Wrapper for adding a ncpd test (python driver)
macro(add_ncpd_test input)
  add_driver_test(${ARGV} EXEC $<TARGET_FILE:Exe_ncpd>)
endmacro()

# ------------------------------------------------------------------------------------------------ #
# end add_mcgrid_test.cmake
# ------------------------------------------------------------------------------------------------ #
