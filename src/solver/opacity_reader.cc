//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/opacity_reader.cc
 * \author <user>
 * \date   <date>
 * \brief  <start>
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "opacity_reader.hh"

namespace odd_solver {

opacity_reader::opacity_reader(const std::string &ipcressfile)
    : my_opacity_file(ipcressfile) {}

void opacity_reader::print_file_info() const { my_opacity_file.printSummary(); }

void opacity_reader::print_mat_info(const size_t mat_id) const {
  my_opacity_file.printSummary(mat_id);
}

} // namespace odd_sovler

//------------------------------------------------------------------------------------------------//
// end of opacity_reader.cc
//------------------------------------------------------------------------------------------------//

