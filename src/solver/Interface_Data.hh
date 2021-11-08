//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/Interface_Data.hh
 * \author Mathew Cleveland
 * \brief  Define the Solver Interface data class that holds the fundamental interface data in
 * convenient c++ data types (ie vectors rather then pointers)
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#ifndef solver_Interface_Data_hh
#define solver_Interface_Data_hh
#include <vector>

namespace odd_solver {

//================================================================================================//
/*!
 * \class Interface_Data
 * \brief
 *
 * This holds the fundamental interface data in simple C++ data types
 *
 */
//================================================================================================//

class Interface_Data {
public:
  //! Default constructors.
  Interface_Data() : number_of_cells(0), number_of_mats(0){};

  bool valid();

  // DATA
  // Data sizes
  size_t number_of_cells;
  size_t number_of_mats;
  // Raw Material data arrays
  std::vector<size_t> problem_matids;
  std::vector<size_t> number_of_cell_mats;
  std::vector<std::vector<size_t>> cell_mats;
  std::vector<std::vector<double>> cell_mat_vol_frac;
  std::vector<std::vector<double>> cell_mat_temperature;
  std::vector<std::vector<double>> cell_mat_density;
};

} // namespace odd_solver

#endif // solver_Interface_Data_hh

//------------------------------------------------------------------------------------------------//
// end of solver/Interface_Data.hh
//------------------------------------------------------------------------------------------------//
