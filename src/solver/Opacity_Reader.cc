//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/Opacity_Reader.cc
 * \author <user>
 * \date   <date>
 * \brief  <start>
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Opacity_Reader.hh"
#include "cdi/OpacityCommon.hh"

namespace odd_solver {

//================================================================================================//
/*!
 * \brief
 *
 * Default constructor to open and read all opacities in an ipcressfile
 *
 *  \param[in] ipcressfile opacity file name
 *
 */
//================================================================================================//
Opacity_Reader::Opacity_Reader(const std::string &ipcressfile)
    : my_opacity_file(new rtt_cdi_ipcress::IpcressFile(ipcressfile)) {
  for (auto &id : my_opacity_file->getMatIDs()) {
    mat_rosseland_abs_models.push_back(std::make_unique<rtt_cdi_ipcress::IpcressGrayOpacity>(
        my_opacity_file, id, rtt_cdi::Model::ROSSELAND, rtt_cdi::Reaction::ABSORPTION));
    mat_planck_abs_models.push_back(std::make_unique<rtt_cdi_ipcress::IpcressGrayOpacity>(
        my_opacity_file, id, rtt_cdi::Model::PLANCK, rtt_cdi::Reaction::ABSORPTION));
    mat_rosseland_total_models.push_back(std::make_unique<rtt_cdi_ipcress::IpcressGrayOpacity>(
        my_opacity_file, id, rtt_cdi::Model::ROSSELAND, rtt_cdi::Reaction::TOTAL));
  }
}

//================================================================================================//
/*!
 * \brief
 *
 * Constructor to open an ipcress file and read a provided list of materials opacities.
 *
 *  \param[in] ipcressfile opacity file name
 *  \param[in] matids list of material id's to read from the file
 *
 */
//================================================================================================//
Opacity_Reader::Opacity_Reader(const std::string &ipcressfile, const std::vector<size_t> &matids)
    : my_opacity_file(new rtt_cdi_ipcress::IpcressFile(ipcressfile)) {
  for (auto &id : matids) {
    mat_rosseland_abs_models.push_back(std::make_unique<rtt_cdi_ipcress::IpcressGrayOpacity>(
        my_opacity_file, id, rtt_cdi::Model::ROSSELAND, rtt_cdi::Reaction::ABSORPTION));
    mat_planck_abs_models.push_back(std::make_unique<rtt_cdi_ipcress::IpcressGrayOpacity>(
        my_opacity_file, id, rtt_cdi::Model::PLANCK, rtt_cdi::Reaction::ABSORPTION));
    mat_rosseland_total_models.push_back(std::make_unique<rtt_cdi_ipcress::IpcressGrayOpacity>(
        my_opacity_file, id, rtt_cdi::Model::ROSSELAND, rtt_cdi::Reaction::TOTAL));
  }
}

} // namespace odd_solver

//------------------------------------------------------------------------------------------------//
// end of Opacity_Reader.cc
//------------------------------------------------------------------------------------------------//
