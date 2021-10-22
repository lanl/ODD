#--------------------------------------------*-cmake-*---------------------------------------------#
# file   odd/autodoc/CMakeLists.txt
# author Mathew Cleveland
# brief  Instructions for building Doxygen documentation
# note   Copyright (C) 2013-2021 Triad National Security, LLC., All rights reserved.
# ------------------------------------------------------------------------------------------------ #

# Draco's config/vendor_libraries.cmake is responsible for detecting if doxygen is available.

if(NOT DOXYGEN_FOUND)
  return()
endif()

include(autodoc_macros)
set_autodocdir() # --> AUTODOCDIR

# Source directories that should be examined by doxygen Sets ${DOXYGEN_INPUT}. This value is
# embedded in the doxygen_config file below.
set_doxygen_input() # --> DOXYGEN_INPUT, DOXYGEN_EXAMPLE_PATH


separate_arguments(doxy_dirs UNIX_COMMAND "${DOXYGEN_INPUT} ${DOXYGEN_EXPAMPLE_PATH}")
foreach(dir ${doxy_dirs})
  if(NOT IS_DIRECTORY ${dir})
    message(
      FATAL_ERROR "Directory ${dir} listed as source path for doxygen, but it does not exist.")
  endif()
endforeach()
unset(doxy_dirs)

# Tell doxygen where image files are located so they can be copied to the output directory.
set_doxygen_image_path() # --> DOXYGEN_IMAGE_PATH

# Doxygen only allows 32 threads max
set_doxygen_dot_num_threads() # --> DOXYGEN_DOT_NUM_THREADS, DOXYGEN_HTML_DYNAMIC_MENUS,
                              # DOXYGEN_WARN_AS_ERROR

set(DOXYGEN_HTML_OUTPUT odd)
set(project_brief "Orthogonal Diffusion Discretization (odd) library.")
string(REPLACE " " "___" project_brief ${project_brief})
string(REPLACE " " "___" DOXYGEN_HTML_DYNAMIC_MENUS "${DOXYGEN_HTML_DYNAMIC_MENUS}")

set_doxygen_tagfiles() # --> TAGFILES
if(TAGFILES MATCHES NOTFOUND)
  unset(TAGFILES)
else()
  set(TAGFILES "-DTAGFILES=\"${TAGFILES}\"")
endif()

# doxygen_config.in needs the following symbols:
add_custom_command(
  OUTPUT "${PROJECT_BINARY_DIR}/autodoc/doxygen_config"
  COMMAND
    "${CMAKE_COMMAND}" -DINFILE="${draco_DIR}/doxygen_config.in"
    -DOUTFILE="${PROJECT_BINARY_DIR}/autodoc/doxygen_config" -DPROJECT_NAME=${PROJECT_NAME}
    -Dpackage=${PROJECT_NAME} -DPROJECT_NUMBER="${Odd_VERSION_FULL}"
    -Dnumber="${Odd_VERSION_FULL}" -DOUTPUT_DIRECTORY="${DOXYGEN_OUTPUT_DIR}"
    -Ddoxygen_output_top="${DOXYGEN_OUTPUT_DIR}" -DINPUT=${DOXYGEN_INPUT}
    -Ddoxygen_input=${DOXYGEN_INPUT} -DDOXYGEN_EXAMPLE_PATH=${DOXYGEN_EXAMPLE_PATH}
    -Ddoxygen_examples=${DOXYGEN_EXAMPLES} -DHTML_OUTPUT=${DOXYGEN_HTML_OUTPUT}
    -Ddoxygen_html_output=${DOXYGEN_HTML_OUTPUT} -DDOTFILE_DIRS="" -Ddotpath=""
    -Dproject_brief=\"${project_brief}\" -Dlatex_yes_no=NO
    -Ddoxygen_latex_output=${DOXYGEN_HTML_OUTPUT}/../latex
    -Dhtmlheader="${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/header.html"
    -Dhtmlfooter="${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/footer.html"
    -Dproject_logo="${PROJECT_SOURCE_DIR}/autodoc/html/montecarlo.jpg"
    -Dnumthreads=${DOXYGEN_DOT_NUM_THREADS} -DDOXYGEN_IMAGE_PATH=${DOXYGEN_IMAGE_PATH} "${TAGFILES}"
    -DDOXYGEN_HTML_DYNAMIC_MENUS="${DOXYGEN_HTML_DYNAMIC_MENUS}"
    -DDOXYGEN_WARN_AS_ERROR="${DOXYGEN_WARN_AS_ERROR}"
    -DDOXYGEN_INCLUDE_PATH="${DOXYGEN_INCLUDE_PATH}"
    -DDOXYGEN_ENABLED_SECTIONS="${DOXYGEN_ENABLED_SECTIONS}" -P
    "${draco_DIR}/configureFileOnMake.cmake"
  DEPENDS "${draco_DIR}/doxygen_config.in" "${PROJECT_SOURCE_DIR}/autodoc/CMakeLists.txt"
          "${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/header.html"
          "${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/footer.html"
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/autodoc
  COMMENT "Generating doxygen_config file.")

# Create header.html, footer.html and doxygen.css.
doxygen_provide_support_files()

add_custom_command(
  OUTPUT "${PROJECT_BINARY_DIR}/autodoc/mainpage.dcc"
  COMMAND
    "${CMAKE_COMMAND}" -DINFILE="${PROJECT_SOURCE_DIR}/autodoc/mainpage.dcc.in"
    -DOUTFILE="${PROJECT_BINARY_DIR}/autodoc/mainpage.dcc" -DCOMP_LINKS=${COMP_LINKS}
    -DPACKAGE_LINKS=${PACKAGE_LINKS} -DOdd_VERSION_FULL=${Odd_VERSION_FULL}
    -DPROJECT_SOURCE_DIR="${PROJECT_SOURCE_DIR}" -P
    "${PROJECT_SOURCE_DIR}/autodoc/generate_mainpage_dcc.cmake"
  DEPENDS "${PROJECT_SOURCE_DIR}/autodoc/mainpage.dcc.in" Exe_imc_info
  COMMENT "Generating mainpage.dcc.")

# The actual doxygen command that will be run
add_custom_command(
  OUTPUT ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/always_out_of_date.html
  BYPRODUCTS ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/index.html
  COMMAND "${DOXYGEN_EXECUTABLE}" ${PROJECT_BINARY_DIR}/autodoc/doxygen_config
  DEPENDS ${PROJECT_BINARY_DIR}/autodoc/doxygen_config ${PROJECT_BINARY_DIR}/autodoc/mainpage.dcc
          ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/doxygen.css
  WORKING_DIRECTORY ${PROJECT_BINARY_DIR}/autodoc
  COMMENT "Generating HTML developer documentation with doxygen.")

# Target for building doxygen files.
add_custom_target(
  autodoc
  DEPENDS ${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/always_out_of_date.html
  COMMENT
    "Building Doxygen based autodocs (HTML)...

*** Autodoc HTML pages are at
    file://${DOXYGEN_OUTPUT_DIR}/${DOXYGEN_HTML_OUTPUT}/index.html

    ")

# This target should be defined by src/micgrid/python/doc/CMakeLists.txt, but only if sphinx-build
# and sphinx-apidoc binaries are discovered.
if(TARGET mcgrid_parser_html_docs)
  add_dependencies(autodoc mcgrid_parser_html_docs)
endif()

# To Do: 1.  How to 'make clean'?

# ------------------------------------------------------------------------------------------------ #
# End autodoc/CMakeLists.txt
# ------------------------------------------------------------------------------------------------ #
