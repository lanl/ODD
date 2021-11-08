//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   api/Arguments.hh
 * \author Mathew Cleveland
 * \brief  API argument definitions 
 * \note   Copyright (C) 2021 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef api_arguments_hh
#define api_arguments_hh

#include <string>

extern "C" {

struct Control_Data {
  std::string opacity_file;
  // Constructor to initialize data
  Control_Data();
  void check_arguments();
};

struct Zonal_Data {
  // Mesh data
  size_t number_of_cells;
  size_t dimensions;
  double dx;
  double dy;
  double dz;

  // Material data
  size_t number_of_mats;
  size_t *problem_matids;
  size_t *number_of_cell_mats;
  size_t *cell_mats;
  double *cell_mat_vol_frac;
  double *cell_mat_temperature;
  double *cell_mat_density;

  // constructor to initialize data
  Zonal_Data();
  void check_arguments();
};

struct Output_Data {
  // grey_opacity_data
  double *ave_opacity_data;

  // constructor to initialize data
  Output_Data();
  void check_arguments();
};

struct Arguments {
  Control_Data control_data;
  Zonal_Data zonal_data;
  Output_Data output_data;
  // constructor to initialize all data fields
  Arguments();
};

} // extern "C"

#endif // api_arguments_hh

//------------------------------------------------------------------------------------------------//
// end of api/arguments.hh
//------------------------------------------------------------------------------------------------//
