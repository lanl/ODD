//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   driver/odd.cc
 * \author Mathew Cleveland
 * \date   October 21st 2021
 * \brief  Example driver
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "api/Arguments.hh"
#include "api/Function_Interface.hh"
#include "driver/Odd_Functions.hh"
#include "c4/global.hh"
#include "ds++/dbc.hh"
#include <iostream>
#include <string>
#include <vector>

int main(int argc, char *argv[]) {

  // Initialize MPI
  rtt_c4::initialize(argc, argv);

  if (argc < 2) {
    std::cout << "Error: Missing Arguments" << std::endl;
    exit(EXIT_FAILURE);
  }

  std::vector<std::string> argv_strings;
  argv_strings.assign(argv + 1, argv + argc);

  // setup data arrays
  Arguments arg;
  odd_driver::Odd_Driver_Data odd_data;

  // Populate data
  odd_driver::build_arguments_from_cmd(argv_strings, arg, odd_data);

  double t = arg.control_data.dt;
  for (size_t cycle = 0; cycle < odd_data.n_cycles; cycle++) {
    const bool print_cycle = ((cycle + 1) % odd_data.print_frequency == 0);

    arg.control_data.print = print_cycle ? 1 : 0;
    // Call the solver on the fake arguments list
    Odd_Diffusion_Solve(arg);

    if (print_cycle) {
      // print out the opacity data
      std::cout << "####################################\n";
      std::cout << "    t = " << t << "\n";
      std::cout << "    cycle = " << cycle + 1 << "\n";
      std::cout << "    dt = " << arg.control_data.dt << "\n";
      std::cout << "####################################\n";
    }

    // Calculate the energy update
    energy_update(arg, odd_data, print_cycle);

    if (print_cycle) {
      // print out the final material and cell properties
      size_t mat_index = 0;
      for (size_t i = 0; i < arg.zonal_data.number_of_local_cells; i++) {
        std::cout << "Cell = " << i;
        std::cout << " loc = (" << arg.zonal_data.cell_position[i * 3];
        for (size_t d = 1; d < arg.zonal_data.dimensions; d++)
          std::cout << ", " << arg.zonal_data.cell_position[i * 3 + d];
        std::cout << ")";
        std::cout << " erad[" << i << "] = " << arg.output_data.cell_erad[i];
        std::cout << " Trad[" << i << "] = " << arg.output_data.cell_Trad[i];
        for (size_t m = 0; m < arg.zonal_data.number_of_cell_mats[i]; m++, mat_index++) {
          std::cout << " cell_mat_delta_e[" << i << "][" << m
                    << "] = " << arg.output_data.cell_mat_delta_e[mat_index];
          std::cout << " cell_mat_T[" << i << "][" << m
                    << "] = " << arg.zonal_data.cell_mat_temperature[mat_index];
        }
        std::cout << "\n";
      }
      std::cout << "####################################\n" << std::endl;
    }
    for (size_t i = 0; i < arg.zonal_data.number_of_local_cells; i++)
      arg.zonal_data.cell_erad[i] = arg.output_data.cell_erad[i];
    t += arg.control_data.dt;
  }
}

//------------------------------------------------------------------------------------------------//
// end of driver.cc
//------------------------------------------------------------------------------------------------//
