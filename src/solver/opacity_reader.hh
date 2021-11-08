//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/opacity_reader.hh
 * \author Mathew Cleveland
 * \brief  Define class opacity_reader
 * \note   Copyright (C) 2021 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef solver_opacity_reader_hh
#define solver_opacity_reader_hh

#include "cdi_ipcress/IpcressFile.hh"
#include "cdi_ipcress/IpcressGrayOpacity.hh"

namespace odd_solver {

//================================================================================================//
/*!
 * \class opacity_reader
 * \brief
 *
 * Example of how to import draco functions
 *
 * \sa opacity_reader.cc for detailed descriptions.
 *
 */
/*!
 * \example solver/test/tstopacity_reader.cc
 *
 * Test of opacity_reader.
 */
//================================================================================================//

class opacity_reader {
public:
  // NESTED CLASSES AND TYPEDEFS

  // CREATORS

  //! Default constructors.
  opacity_reader(const std::string &ipcressFile);

  //! Load up the requested matids
  opacity_reader(const std::string &ipcressFile, const std::vector<size_t> &matids);

  // SERVICES
  void print_available_mats() const {
    std::cout << "Available matids " << std::endl;
    std::cout << "    ";
    for (auto &id : my_opacity_file->getMatIDs())
      std::cout << std::to_string(id) << " ";
    std::cout << std::endl;
  }

  void print_available_data(const size_t mat_id) const {
    std::cout << "Available ipcress fields " << std::endl;
    std::cout << "    ";
    for (auto &string : my_opacity_file->listDataFieldNames(mat_id))
      std::cout << string << " ";
    std::cout << std::endl;
  }

private:
  // DATA
  std::shared_ptr<rtt_cdi_ipcress::IpcressFile> my_opacity_file;

public:
  std::vector<std::unique_ptr<rtt_cdi_ipcress::IpcressGrayOpacity>> mat_rosseland_abs_models;
};

} // namespace odd_solver

#endif // solver_opacity_reader_hh

//------------------------------------------------------------------------------------------------//
// end of solver/opacity_reader.hh
//------------------------------------------------------------------------------------------------//
