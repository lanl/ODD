//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   driver/driver.cc
 * \author Mathew Cleveland
 * \date   October 21st 2021
 * \brief  Example driver
 * \note   Copyright (C) 2021 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "api/Arguments.hh"
#include "api/Function_Interface.hh"
#include <iostream>
#include <vector>

int main(int argc, char *argv[]) {
  Arguments arg;

  // fill control data with "valid junk" and check arguments
  arg.control_data.opacity_file = argv[1];
  arg.control_data.check_arguments();

  // fill the zonal data with "valid junk" and check arguments
  arg.zonal_data.number_of_cells = 2;
  arg.zonal_data.dimensions = 1;
  arg.zonal_data.dx = 31.0;
  // global material data
  std::vector<size_t> matids = {10001, 10002};
  arg.zonal_data.number_of_mats = 2;
  arg.zonal_data.problem_matids = &matids[0];
  // cell wise material data
  // 1 material in cell 1 and 2 materials in cell 2
  std::vector<size_t> cell_number_of_mats{1, 2};
  arg.zonal_data.number_of_cell_mats = &cell_number_of_mats[0];
  // cell 1 (mat 1) cell 2 (mat 1 and mat 2)
  std::vector<size_t> cell_mats{1, 1, 2};
  // cell 1 (1.0) cell 2 mat 1 (0.5) and mat 2 (0.5)
  std::vector<double> cell_mat_vol_frac{1.0, 0.5, 0.5};
  // cell 1 (1.0) cell 2 mat 1 (10.0) and mat 2 (1.0)
  std::vector<double> cell_mat_temperature{1.0, 10.0, 1.0};
  // cell 1 (10) cell 2 mat 1 (1.0) and mat 2 (10.0)
  std::vector<double> cell_mat_density{10.0, 1.0, 10.0};
  arg.zonal_data.cell_mats = &cell_mats[0];
  arg.zonal_data.cell_mat_vol_frac = &cell_mat_vol_frac[0];
  arg.zonal_data.cell_mat_temperature = &cell_mat_temperature[0];
  arg.zonal_data.cell_mat_density = &cell_mat_density[0];

  // setup the opacity data field to be filled in by the solver
  std::vector<double> opacity_data(arg.zonal_data.number_of_cells, 0.0);
  arg.output_data.ave_opacity_data = &opacity_data[0];

  // Call the solver on the fake arguments list
  std::cout << "Call solver" << std::endl;
  Odd_Diffusion_Solve(arg);

  // print out the opacity data
  auto optr = arg.output_data.ave_opacity_data;
  for (size_t i = 0; i < arg.zonal_data.number_of_cells; i++, optr++)
    std::cout << "Opacity_value[" << i << "] = " << *optr << std::endl;
}

//------------------------------------------------------------------------------------------------//
// end of <basename>.cc
//------------------------------------------------------------------------------------------------//
