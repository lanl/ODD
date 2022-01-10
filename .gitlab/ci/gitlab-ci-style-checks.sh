#!/bin/bash -l
#--------------------------------------------------------------------------------------------------#
# File  : ./.gitlab/ci/gitlab-style-checks.sh
# Date  : Tuesday, Jun 02, 2020, 12:28 pm
# Author: Kelly Thompson <kgt@lanl.gov>
# Note  : Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
#--------------------------------------------------------------------------------------------------#

# preliminaries and environment
set -e
#shellcheck source=.gitlab/ci/common.sh
source "${ODD_SOURCE_DIR}/.gitlab/ci/common.sh"
#shellcheck source=.gitlab/ci/environments.sh
source "${ODD_SOURCE_DIR}/.gitlab/ci/environments.sh"

echo -e "\n========== printenv ==========\n"
if ! [[ "${SLURM_NODELINST:-notset}" == "notset" ]]; then
  echo "SLURM_NODELIST = ${SLURM_NODELIST}"
fi
echo "HOSTNAME       = ${HOSTNAME}"
# printenv
# run "module avail"
printenv | grep "CI_"
echo " "

#--------------------------------------------------------------------------------------------------#
# Style check
#--------------------------------------------------------------------------------------------------#

run "cd $ODD_SOURCE_DIR"
run "git fetch --all"
run "git branch -a"
if [[ $(git branch | grep -c master) == 0 ]]; then
  run "git branch master origin/master"
fi
run "${ODD_SOURCE_DIR}/.gitlab/ci/check_style.sh" -t

echo -e "\n======== end .gitlab-ci-run-tests.sh ==========\n"

#--------------------------------------------------------------------------------------------------#
# End .gitlab/ci/gitlab-ci-run-tests.sh
#--------------------------------------------------------------------------------------------------#
