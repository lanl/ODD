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
#include "c4/opstream.hh"
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

  double t = 0.0;
  for (size_t cycle = 0; cycle < odd_data.n_cycles; cycle++) {
    const bool print_cycle = ((cycle + 1) % odd_data.print_frequency == 0);

    // calculate volume source if specified
    if (odd_data.vol_source_strength > 0.0)
      odd_driver::update_source(arg, odd_data, t);

    arg.control_data.print = print_cycle ? 1 : 0;
    // Call the solver on the fake arguments list
    Odd_Diffusion_Solve(arg);

    rtt_c4::opstream p_out;
    p_out.precision(4);
    p_out.setf(std::ios::scientific, std::ios::floatfield);
    if (print_cycle && rtt_c4::node() == 0) {
      // print out the opacity data
      p_out << "############ CYCLE BEGIN ############\n";
      p_out << "    t = " << t + arg.control_data.dt << "\n";
      p_out << "    cycle = " << cycle + 1 << "\n";
      p_out << "    dt = " << arg.control_data.dt << "\n";
      p_out << "####################################\n";
    }
    p_out.send();
    // Calculate the energy update
    energy_update(arg, odd_data, print_cycle);

    if ((odd_data.domain_decomposed && print_cycle) || (print_cycle && rtt_c4::node() == 0)) {

      // print out the final material and cell properties
      size_t mat_index = 0;
      size_t group_index = 0;
      for (size_t i = 0; i < arg.zonal_data.number_of_local_cells; i++) {
        const size_t cell_id = arg.zonal_data.cell_global_id[i];
        p_out << "Cell = " << cell_id;
        p_out << "; loc = (" << arg.zonal_data.cell_position[i * 3];
        for (size_t d = 1; d < arg.zonal_data.dimensions; d++)
          p_out << ", " << arg.zonal_data.cell_position[i * 3 + d];
        p_out << ")";
        p_out << "; erad[" << cell_id << "] = " << arg.output_data.cell_erad[i];
        p_out << "; Trad[" << cell_id << "] = " << arg.output_data.cell_Trad[i];
        for (size_t m = 0; m < arg.zonal_data.number_of_cell_mats[i]; m++, mat_index++) {
          p_out << "; cell_mats[" << cell_id << "][" << m
                << "] = " << arg.zonal_data.cell_mats[mat_index];
          p_out << "; cell_mat_delta_e[" << cell_id << "][" << m
                << "] = " << arg.output_data.cell_mat_delta_e[mat_index];
          p_out << "; cell_mat_T[" << cell_id << "][" << m
                << "] = " << arg.zonal_data.cell_mat_temperature[mat_index];
        }
        for (size_t g = 0; g < arg.control_data.ngroups; g++, group_index++) {
          p_out << "; mg_erad[" << cell_id << "][" << g
                << "] = " << arg.zonal_data.cell_mg_erad[group_index];
        }
        if (odd_data.domain_decomposed)
          p_out << "; Rank = " << rtt_c4::node();
        p_out << "\n";
      }
      if (odd_data.domain_decomposed)
        p_out.send();
      if (rtt_c4::node() == 0)
        p_out << "############  CYCLE END  ############\n";
    }
    p_out.send();
    t += arg.control_data.dt;
    arg.control_data.dt *= odd_data.dt_ramp;
  }
}

//------------------------------------------------------------------------------------------------//
// end of driver.cc
//------------------------------------------------------------------------------------------------//
