#!/bin/bash -l
#--------------------------------------------------------------------------------------------------#
# File  : ./.gitlab/ci/gitlab-ci-run-tests.sh
# Date  : Tuesday, Jun 02, 2020, 12:28 pm
# Author: Kelly Thompson <kgt@lanl.gov>
# Note  : Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

# ODD protect for now with ccsrad group
umask 0007
umask

# preliminaries and environment
set -e
# shellcheck source=.gitlab/ci/common.sh
source "${ODD_SOURCE_DIR}/.gitlab/ci/common.sh"
# shellcheck source=.gitlab/ci/environments.sh
source "${ODD_SOURCE_DIR}/.gitlab/ci/environments.sh"

if [[ $(which lscpu | grep -c lscpu) -gt 0 ]]; then
  NPROC=$(grep -c processor < /proc/cpuinfo)
fi

[[ -z "${CTEST_NPROC}" ]] && CTEST_NPROC=$NPROC || echo "limiting CTEST_NPROC = $CTEST_NPROC"
[[ -z "${MAXLOAD}" ]] && MAXLOAD=$NPROC || echo "limiting MAXLOAD = $MAXLOAD"
[[ "${DEPLOY}" == "TRUE" ]] && EXTRA_CMAKE_ARGS="-DBUILD_TESTING=NO ${EXTRA_CMAKE_ARGS}"
if [[ -n "${ODD_INSTALL_DIR}" ]]; then
  if ! [[ -d "${ODD_INSTALL_DIR}" ]]; then
    run "mkdir -p ${ODD_INSTALL_DIR}"
    run "chgrp ccsrad ${ODD_INSTALL_DIR}"
    run "chmod g+rwX,o-rX ${ODD_INSTALL_DIR}"
    run "chmod g+s ${ODD_INSTALL_DIR}"
  fi
   EXTRA_CMAKE_ARGS="-DCMAKE_INSTALL_PREFIX=${ODD_INSTALL_DIR} ${EXTRA_CMAKE_ARGS}"
fi

# We expect draco to be found at DRACO_INSTALL_DIR
[[ -z "${DRACO_INSTALL_DIR}" ]] && die "DRACO_INSTALL_DIR is not set."
[[ ! -d "${DRACO_INSTALL_DIR}" ]] && die "DRACO_INSTALL_DIR=${DRACO_INSTALL_DIR} not found."
CMAKE_PREFIX_PATH="${CMAKE_PREFIX_PATH:+"$CMAKE_PREFIX_PATH:"}${DRACO_INSTALL_DIR}"
export CMAKE_PREFIX_PATH

#if [[ -d /ccs/opt/texlive/2018/texmf-dist/tex/latex/newunicodechar ]]; then
#  TEXINPUTS="${TEXINPUTS:+${TEXINPUTS}:}/ccs/opt/texlive/2018/texmf-dist/tex/latex/newunicodechar"
#fi

# Report number of cores/threads, etc.
# run "lstopo --no-io --no-bridges --no-icaches"

echo -e "\n========== printenv ==========\n"
[[ -z "${SLURM_NODELIST}" ]] || echo "SLURM_NODELIST = ${SLURM_NODELIST}"
echo "HOSTNAME       = ${HOSTNAME}"
echo -e "NPROC       = ${NPROC}\n"
echo -e "CTEST_NPROC = ${CTEST_NPROC}\n"
echo -e "MAXLOAD     = ${MAXLOAD}\n"
echo -e "EXTRA_CMAKE_ARGS = ${EXTRA_CMAKE_ARGS}\n"
echo -e "CMAKE_BUILD_TYPE = ${CMAKE_BUILD_TYPE}\n"
echo -e "EXTRA_CTEST_ARGS = ${EXTRA_CTEST_ARGS}\n"
echo -e "DRACO_INSTALL_DIR = ${DRACO_INSTALL_DIR}\n"
# echo -e "CMAKE_PREFIX_PATH = ${CMAKE_PREFIX_PATH}\n"
if [[ "${AUTODOC}" == "ON" ]]; then
  run "which doxygen"
  run "which latex"
  run "which pdflatex"
  echo "TEXINPUTS = ${TEXINPUTS}"
fi
if [[ "${EXTRA_CTEST_ARGS}" =~ memcheck ]]; then
  run "which valgrind"
fi
printenv >& environment.log
run "pwd"

[[ "${EXTRA_CMAKE_ARGS}" =~ "CODE_COVERAGE" ]] && CODECOV=ON || CODECOV=OFF

#--------------------------------------------------------------------------------------------------#
# Setup compiler flags
#--------------------------------------------------------------------------------------------------#

for i in C CXX Fortran CUDA; do

  # Turn build warnings into fatal errors
  if [[ "${DRACO_ENV}" =~ "xl16" ]]; then
    eval export ${i}_FLAGS+=\" -qhalt=w\"
  else
    if [[ "${i}" == "Fortran" ]] && [[ "${FC}" =~ "ifort" ]]; then
      # do not add the flag
      echo "Skip adding -Werror to ${i}_FLAGS (FC = ${FC})"
    elif [[ "${i}" == "CUDA" ]]; then
      eval export ${i}_FLAGS+=\" -Werror all-warnings\"
    else
      # Works for gcc, llvm, and icpc (not ifort).
      eval export ${i}_FLAGS+=\" -Werror\"
    fi
  fi

  # Enable coverage for Debug builds
  if [[ ${CMAKE_BUILD_TYPE} == Debug ]] && [[ "${CODECOV}" == "ON" ]] ; then
    eval export ${i}_FLAGS+=\" --coverage\"
  fi
done
echo "C_FLAGS       = $C_FLAGS"
echo "CXX_FLAGS     = $CXX_FLAGS"
# shellcheck disable=2154
echo "Fortran_FLAGS = $Fortran_FLAGS"
echo "CUDA_FLAGS = $CUDA_FLAGS"

# Tweak MPI settings
# export OMPI_MCA_btl=self,sm
export OMPI_MCA_btl=^openib

#--------------------------------------------------------------------------------------------------#
# Odd
#--------------------------------------------------------------------------------------------------#

echo -e "\n========== Odd ==========\n"
if ! [[ -d ${ODD_BINARY_DIR} ]]; then run "mkdir -p ${ODD_BINARY_DIR}"; fi
run "cd ${ODD_BINARY_DIR}"
echo -e "\n========"
config_cmd="cmake ${EXTRA_CMAKE_ARGS} -DCMAKE_BUILD_TYPE=${CMAKE_BUILD_TYPE}"
config_cmd+=" ${ODD_SOURCE_DIR}"
echo "==> $config_cmd"
eval "$config_cmd"
echo " "
echo -e "\n========\n"
if [[ "${AUTODOC}" == "ON" ]]; then
  run "make -j -l $MAXLOAD autodoc"
  # Run one test so that the CDash report isn't empty.
  run "ctest -R odd_release_tstOddRelease"
else
  run "make -j -l $MAXLOAD" #  VERBOSE=1
  echo -e "\n========\n"
  if [[ "${TEST_EXCLUSIONS:-no}" != "no"  ]]; then
    EXTRA_CTEST_ARGS+=" -E ${TEST_EXCLUSIONS}"
  fi
  if [[ "${DEPLOY}" == "TRUE" ]]; then
    run "make -j -l $MAXLOAD install"
  else
    test_cmd="ctest ${EXTRA_CTEST_ARGS} -j ${CTEST_NPROC} --test-load ${MAXLOAD}"
    test_cmd+=" --output-on-failure --stop-on-failure"
    echo "==> $test_cmd"
    eval "$test_cmd"
  fi
fi

# Generate a coverage report (text and html)
if [[ "${CODECOV}" == "ON" && "${DEPLOY}" != "TRUE" ]]; then
  run "make covrep"
fi

# Double check permissions for installed (deployed) files
if [[ "${DEPLOY}" == "TRUE" ]]; then
  # run "find ${ODD_INSTALL_DIR} ! -perm -g+rw -exec ls -aFl {} \;"
  # badpermfiles=$(find ${ODD_INSTALL_DIR} ! -perm -g+rw -exec ls  {} \;)
  run "chgrp -R ccsrad ${ODD_INSTALL_DIR}"
  run "chmod -R g+rwX,o-rX ${ODD_INSTALL_DIR}"
fi

if [[ ${CTEST_MODE} == "Nightly" ]]; then
  echo -e "\n========\n"
  echo "Nighlty regressions should run a ctest in script mode and post to CDash"
fi

echo -e "\n======== end .gitlab-ci-run-tests.sh ==========\n"

#--------------------------------------------------------------------------------------------------#
# End .gitlab/ci/gitlab-ci-run-tests-ccsnet.sh
#--------------------------------------------------------------------------------------------------#
