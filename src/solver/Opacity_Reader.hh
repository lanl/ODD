//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/Opacity_Reader.hh
 * \author Mathew Cleveland
 * \brief  Define class Opacity_Reader
 * \note   Copyright (C) 2021 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef solver_Opacity_Reader_hh
#define solver_Opacity_Reader_hh

#include "cdi_ipcress/IpcressFile.hh"
#include "cdi_ipcress/IpcressGrayOpacity.hh"

namespace odd_solver {

//================================================================================================//
/*!
 * \class Opacity_Reader
 * \brief
 *
 * Example of how to import draco functions
 *
 * \sa Opacity_Reader.cc for detailed descriptions.
 *
 */
/*!
 * \example solver/test/tstOpacity_Reader.cc
 *
 * Test of Opacity_Reader.
 */
//================================================================================================//

class Opacity_Reader {
public:
  // NESTED CLASSES AND TYPEDEFS

  // CREATORS

  //! Default constructors.
  Opacity_Reader(const std::string &ipcressFile);

  //! Load up the requested matids
  Opacity_Reader(const std::string &ipcressFile, const std::vector<size_t> &matids);

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
  std::vector<std::unique_ptr<rtt_cdi_ipcress::IpcressGrayOpacity>> mat_planck_abs_models;
  std::vector<std::unique_ptr<rtt_cdi_ipcress::IpcressGrayOpacity>> mat_rosseland_total_models;
};

} // namespace odd_solver

#endif // solver_Opacity_Reader_hh

//------------------------------------------------------------------------------------------------//
// end of solver/Opacity_Reader.hh
//------------------------------------------------------------------------------------------------//
