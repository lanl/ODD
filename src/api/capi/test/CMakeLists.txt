# -------------------------------------------*-cmake-*-------------------------------------------- #
# file   driver/CMakeLists.txt
# author Mathew Cleveland
# date   October 2021
# brief  Generate build project files for mc
# note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #
project(capi_test C)

add_component_executable(
  TARGET Exe_example
  TARGET_DEPS "Lib_api"
  SOURCES "example.c"
  EXPORT_NAME "odd-targets")

# ------------------------------------------------------------------------------------------------ #
# Installation instructions
# ------------------------------------------------------------------------------------------------ #
install(TARGETS Exe_example
    EXPORT odd-targets
    DESTINATION ${DBSCFGDIR}capi/test)

file(GLOB ipcress_files ${PROJECT_SOURCE_DIR}/*.ipcress ${PROJECT_SOURCE_DIR}/*.gold)
provide_aux_files(FILES "${ipcress_files}" FOLDER capi_test)


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
