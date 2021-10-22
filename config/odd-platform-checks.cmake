# -------------------------------*-cmake-*-------------------------------------------------------- #
# file   : odd-platform-checks.cmake
# brief  : Platform Checks for the Odd project
# note   : Copyright (C) 2018-2021 Triad National Security, LLC., All rights reserved
# ------------------------------------------------------------------------------------------------ #

# Perform extra platform checks that Odd requries that Draco doesn't require.
macro(odd_platform_checks)

if(NOT DEFINED ODD_PLATFORM_CHECKS_DONE)
    set(ODD_PLATFORM_CHECKS_DONE
        TRUE
        CACHE BOOL "Have we completed the platform checks for Odd?")

    message("
Platform Checks (Odd)...
")

    include(CheckIncludeFiles)
    check_include_files(emmintrin.h HAVE_SSE2)

  endif()

endmacro()

# ------------------------------------------------------------------------------------------------ #
# End config/odd-platform-checks.cmake
# ------------------------------------------------------------------------------------------------ #
