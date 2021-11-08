//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/opacity_reader.cc
 * \author <user>
 * \date   <date>
 * \brief  <start>
 * \note   Copyright (C) 2021 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "opacity_reader.hh"
#include "cdi/OpacityCommon.hh"

namespace odd_solver {

opacity_reader::opacity_reader(const std::string &ipcressfile)
    : my_opacity_file(new rtt_cdi_ipcress::IpcressFile(ipcressfile)) {
  for (auto &id : my_opacity_file->getMatIDs())
    mat_rosseland_abs_models.push_back(std::make_unique<rtt_cdi_ipcress::IpcressGrayOpacity>(
        my_opacity_file, id, rtt_cdi::Model::ROSSELAND, rtt_cdi::Reaction::ABSORPTION));
}

opacity_reader::opacity_reader(const std::string &ipcressfile, const std::vector<size_t> &matids)
    : my_opacity_file(new rtt_cdi_ipcress::IpcressFile(ipcressfile)) {
  for (auto &id : matids)
    mat_rosseland_abs_models.push_back(std::make_unique<rtt_cdi_ipcress::IpcressGrayOpacity>(
        my_opacity_file, id, rtt_cdi::Model::ROSSELAND, rtt_cdi::Reaction::ABSORPTION));
}

} // namespace odd_solver

//------------------------------------------------------------------------------------------------//
// end of opacity_reader.cc
//------------------------------------------------------------------------------------------------//
