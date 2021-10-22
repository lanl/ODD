//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/opacity_reader.hh
 * \author Mathew Cleveland
 * \brief  Define class opacity_reader
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef solver_opacity_reader_hh
#define solver_opacity_reader_hh

#include "cdi_ipcress/IpcressFile.hh"

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

    // SERVICES
    void print_file_info() const;
    void print_mat_info(const size_t mat_id) const;

  private:
    // DATA
    rtt_cdi_ipcress::IpcressFile my_opacity_file;
  };

} // end namespace <namespace>

#endif // solver_opacity_reader_hh

//------------------------------------------------------------------------------------------------//
// end of solver/opacity_reader.hh
//------------------------------------------------------------------------------------------------//

