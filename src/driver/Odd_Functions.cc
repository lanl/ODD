//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file  driver/Odd_Functions.cc
 * \author Mathew Cleveland
 * \date   January 21 2022
 * \brief  Basica Interface Function to setup ODD problems
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Odd_Functions.hh"
#include "solver/Constants.hh"
#include "ds++/dbc.hh"
#include <iostream>
#include <string>

namespace odd_driver {

void build_arguments_from_cmd(const std::vector<std::string> argv, Arguments &args,
                              Odd_Driver_Data &odd_data) {

  for (size_t i = 0; i < argv.size(); i++) {
    std::string str = argv[i];
    // Parse control data
    if (str == "-if" || str == "-ipcress_file") {
      odd_data.opacity_file = argv[i + 1];
      args.control_data.opacity_file = &odd_data.opacity_file[0];
    }
    if (str == "-nc" || str == "-n_cycles")
      odd_data.n_cycles = std::stoi(argv[i + 1]);
    if (str == "-dt")
      args.control_data.dt = std::stod(argv[i + 1]);
    if (str == "-mi" || str == "-max_iter")
      args.control_data.max_iter = std::stoi(argv[i + 1]);
    if (str == "-mt" || str == "-min_tol")
      args.control_data.min_tol = std::stod(argv[i + 1]);
    // Mesh data
    if (str == "-cs" || str == "-coord_sys")
      args.zonal_data.coord_sys = std::stoi(argv[i + 1]);
    if (str == "-d" || str == "-dimensions")
      args.zonal_data.dimensions = std::stoi(argv[i + 1]);
    if (str == "-ms" || str == "-mesh_size") {
      odd_data.mesh_size[0] = std::stod(argv[i + 1]);
      odd_data.mesh_size[1] = std::stod(argv[i + 2]);
      odd_data.mesh_size[2] = std::stod(argv[i + 3]);
    }
    if (str == "-mnc" || str == "-mesh_n_cells") {
      odd_data.mesh_n_cells[0] = std::stoi(argv[i + 1]);
      odd_data.mesh_n_cells[1] = std::stoi(argv[i + 2]);
      odd_data.mesh_n_cells[2] = std::stoi(argv[i + 3]);
    }
    // Mat data
    if (str == "-mid" || str == "-matid")
      odd_data.matid = std::stoi(argv[i + 1]);
    if (str == "-tev" || str == "-temperature")
      odd_data.temperature = std::stod(argv[i + 1]);
    if (str == "-rev" || str == "-rad_temperature")
      odd_data.rad_temperature = std::stod(argv[i + 1]);
    if (str == "-rho" || str == "-density")
      odd_data.density = std::stod(argv[i + 1]);
    if (str == "-cv" || str == "-specific_heat")
      odd_data.specific_heat = std::stod(argv[i + 1]);
    if (str == "-bt" || str == "-boundary_temp") {
      odd_data.bnd_temp[0] = std::stod(argv[i + 1]);
      odd_data.bnd_temp[1] = std::stod(argv[i + 2]);
      odd_data.bnd_temp[2] = std::stod(argv[i + 2]);
      odd_data.bnd_temp[3] = std::stod(argv[i + 3]);
      odd_data.bnd_temp[4] = std::stod(argv[i + 4]);
      odd_data.bnd_temp[5] = std::stod(argv[i + 5]);
    }
    if (str == "-rb" || str == "-reflect_bnd") {
      odd_data.reflect_bnd[0] = std::stoi(argv[i + 1]);
      odd_data.reflect_bnd[1] = std::stoi(argv[i + 2]);
      odd_data.reflect_bnd[2] = std::stoi(argv[i + 2]);
      odd_data.reflect_bnd[3] = std::stoi(argv[i + 3]);
      odd_data.reflect_bnd[4] = std::stoi(argv[i + 4]);
      odd_data.reflect_bnd[5] = std::stoi(argv[i + 5]);
    }
  }
  // Assign remaining control data
  args.control_data.bnd_temp = &odd_data.bnd_temp[0];
  args.control_data.reflect_bnd = &odd_data.reflect_bnd[0];

  // Build mesh arguments
  args.zonal_data.domain_decomposed = 0;
  const size_t nfaces_per_cell = args.zonal_data.dimensions * 2;
  size_t ncells = 1;
  for (size_t d = 0; d < args.zonal_data.dimensions; d++) {
    Insist(odd_data.mesh_size[d] > 0, "Included mesh dimension size must be greater then zero");
    ncells *= odd_data.mesh_n_cells[d];
  }
  // allocate/assign mesh data
  args.zonal_data.number_of_local_cells = ncells;
  args.zonal_data.number_of_global_cells = ncells;
  args.zonal_data.number_of_ghost_cells = 0;
  odd_data.cell_position = std::vector<double>(ncells * 3, 0.0);
  args.zonal_data.cell_position = &odd_data.cell_position[0];
  odd_data.cell_size = std::vector<double>(ncells * 3, 0.0);
  args.zonal_data.cell_size = &odd_data.cell_size[0];
  odd_data.cell_global_id = std::vector<size_t>(ncells, 0.0);
  args.zonal_data.cell_global_id = &odd_data.cell_global_id[0];
  odd_data.face_type = std::vector<size_t>(ncells * nfaces_per_cell, 0.0);
  args.zonal_data.face_type = &odd_data.face_type[0];
  odd_data.next_cell_id = std::vector<size_t>(ncells * nfaces_per_cell, 0.0);
  args.zonal_data.next_cell_id = &odd_data.next_cell_id[0];
  // build mesh data
  const std::array<double, 3> delta{
      odd_data.mesh_size[0] / static_cast<double>(odd_data.mesh_n_cells[0]),
      odd_data.mesh_n_cells[1] > 0
          ? odd_data.mesh_size[1] / static_cast<double>(odd_data.mesh_n_cells[1])
          : 0.0,
      odd_data.mesh_n_cells[1] * odd_data.mesh_n_cells[2] > 0
          ? odd_data.mesh_size[2] / static_cast<double>(odd_data.mesh_n_cells[2])
          : 0.0};
  for (size_t cell = 0; cell < ncells; cell++) {
    odd_data.cell_global_id[cell] = cell;
    const std::array<size_t, 3> dim_index{
        cell % odd_data.mesh_n_cells[0],
        odd_data.mesh_n_cells[0] > 0
            ? ((cell / odd_data.mesh_n_cells[0]) %
               (args.zonal_data.dimensions > 2
                    ? odd_data.mesh_n_cells[2]
                    : (odd_data.mesh_n_cells[1] > 0 ? odd_data.mesh_n_cells[1] : 1)))
            : 0,
        odd_data.mesh_n_cells[0] * odd_data.mesh_n_cells[1] > 0
            ? cell / (odd_data.mesh_n_cells[0] * odd_data.mesh_n_cells[1])
            : 0};
    const std::array<size_t, 6> next_cell{
        cell > 0 ? cell - 1 : 0,
        cell + 1,
        cell > odd_data.mesh_n_cells[0] ? cell - odd_data.mesh_n_cells[0] : 0,
        cell + odd_data.mesh_n_cells[0],
        cell > odd_data.mesh_n_cells[0] * odd_data.mesh_n_cells[1]
            ? cell - odd_data.mesh_n_cells[0] * odd_data.mesh_n_cells[1]
            : 0,
        cell + odd_data.mesh_n_cells[0] * odd_data.mesh_n_cells[1]};
    //std::cout << cell << " (";
    for (size_t d = 0; d < args.zonal_data.dimensions; d++) {
      odd_data.cell_size[cell * 3 + d] = delta[d];
      odd_data.cell_position[cell * 3 + d] =
          0.5 * delta[d] + delta[d] * static_cast<double>(dim_index[d]);
      odd_data.face_type[cell * nfaces_per_cell + d * 2] = dim_index[d] == 0 ? 1 : 0;
      odd_data.face_type[cell * nfaces_per_cell + d * 2 + 1] =
          dim_index[d] == (odd_data.mesh_n_cells[d] - 1) ? 1 : 0;
      odd_data.next_cell_id[cell * nfaces_per_cell + d * 2] =
          dim_index[d] == 0 ? ncells : static_cast<size_t>(next_cell[d * 2]);
      odd_data.next_cell_id[cell * nfaces_per_cell + d * 2 + 1] =
          dim_index[d] == (odd_data.mesh_n_cells[d] - 1)
              ? ncells
              : static_cast<size_t>(next_cell[d * 2 + 1]);
      /*
      std::cout << odd_data.cell_size[cell * 3 + d] << " " << odd_data.cell_position[cell * 3 + d]
                << " " << odd_data.face_type[cell * nfaces_per_cell + d * 2] << " "
                << odd_data.face_type[cell * nfaces_per_cell + d * 2 + 1] << " "
                << odd_data.next_cell_id[cell * nfaces_per_cell + d * 2] << " "
                << odd_data.next_cell_id[cell * nfaces_per_cell + d * 2 + 1] << ", ";
                */
    }
    //std::cout << ")" << std::endl;
  }

  // Allocate and assign the material data
  args.zonal_data.number_of_mats = 1;
  odd_data.problem_matids = std::vector<size_t>(1, odd_data.matid);
  args.zonal_data.problem_matids = &odd_data.problem_matids[0];
  odd_data.number_of_cell_mats = std::vector<size_t>(ncells, 1);
  args.zonal_data.number_of_cell_mats = &odd_data.number_of_cell_mats[0];
  odd_data.cell_mats = std::vector<size_t>(ncells, 0);
  args.zonal_data.cell_mats = &odd_data.cell_mats[0];
  odd_data.cell_mat_vol_frac = std::vector<double>(ncells, 1.0);
  args.zonal_data.cell_mat_vol_frac = &odd_data.cell_mat_vol_frac[0];
  odd_data.cell_mat_temperature = std::vector<double>(ncells, odd_data.temperature);
  args.zonal_data.cell_mat_temperature = &odd_data.cell_mat_temperature[0];
  odd_data.cell_mat_density = std::vector<double>(ncells, odd_data.density);
  args.zonal_data.cell_mat_density = &odd_data.cell_mat_density[0];
  odd_data.cell_mat_specific_heat = std::vector<double>(ncells, odd_data.specific_heat);
  args.zonal_data.cell_mat_specific_heat = &odd_data.cell_mat_specific_heat[0];
  odd_data.cell_velocity = std::vector<double>(ncells * 3, 0.0);
  args.zonal_data.cell_velocity = &odd_data.cell_velocity[0];
  odd_data.cell_erad = std::vector<double>(ncells, odd_solver::constants::a *
                                                       std::pow(odd_data.rad_temperature, 4.0));
  args.zonal_data.cell_erad = &odd_data.cell_erad[0];

  // Output data
  odd_data.output_cell_erad = std::vector<double>(
      ncells, odd_solver::constants::a * std::pow(odd_data.rad_temperature, 4.0));
  args.output_data.cell_erad = &odd_data.output_cell_erad[0];
  odd_data.output_cell_Trad = std::vector<double>(
      ncells, odd_solver::constants::a * std::pow(odd_data.rad_temperature, 4.0));
  args.output_data.cell_Trad = &odd_data.output_cell_Trad[0];
  odd_data.output_cell_mat_delta_e = std::vector<double>(ncells, 0.0);
  args.output_data.cell_mat_delta_e = &odd_data.output_cell_mat_delta_e[0];

  size_t mat_index = 0;
  odd_data.total_rad_energy = 0.0;
  odd_data.total_mat_energy = 0.0;
  for (size_t i = 0; i < args.zonal_data.number_of_local_cells; i++) {
    double volume = 1.0;
    for (size_t d = 0; d < args.zonal_data.dimensions; d++)
      volume *= odd_data.cell_size[i * 3 + d];
    odd_data.total_rad_energy += odd_data.cell_erad[i] * volume;
    for (size_t m = 0; m < args.zonal_data.number_of_cell_mats[i]; m++, mat_index++) {
      odd_data.total_mat_energy += odd_data.cell_mat_temperature[mat_index] *
                                   odd_data.cell_mat_specific_heat[mat_index] *
                                   odd_data.cell_mat_vol_frac[mat_index] * volume;
    }
  }
  odd_data.total_energy = odd_data.total_mat_energy + odd_data.total_rad_energy;
}

void energy_update(Arguments &args, Odd_Driver_Data &odd_data) {
  std::cout << " CONSERVATION DATA \n";
  std::cout << "     total rad_energy0 = " << odd_data.total_rad_energy << "\n";
  std::cout << "     total mat_energy0 = " << odd_data.total_mat_energy << "\n";
  std::cout << "     total energy0 = " << odd_data.total_energy << "\n";
  const double total_energy0 = odd_data.total_energy;
  size_t mat_index = 0;
  odd_data.total_rad_energy = 0.0;
  for (size_t i = 0; i < args.zonal_data.number_of_local_cells; i++) {
    double volume = 1.0;
    for (size_t d = 0; d < args.zonal_data.dimensions; d++)
      volume *= odd_data.cell_size[i * 3 + d];
    odd_data.cell_erad[i] = args.output_data.cell_erad[i];
    odd_data.total_rad_energy += odd_data.cell_erad[i] * volume;
    for (size_t m = 0; m < args.zonal_data.number_of_cell_mats[i]; m++, mat_index++) {
      odd_data.total_mat_energy += args.output_data.cell_mat_delta_e[mat_index] * volume;
      // This needs updated for non-constant specific heat
      odd_data.cell_mat_temperature[mat_index] +=
          args.output_data.cell_mat_delta_e[mat_index] /
          (odd_data.cell_mat_specific_heat[mat_index] * odd_data.cell_mat_density[mat_index]);
    }
  }
  odd_data.total_energy = odd_data.total_mat_energy + odd_data.total_rad_energy;

  std::cout << "     total rad_energy = " << odd_data.total_rad_energy << "\n";
  std::cout << "     total mat_energy = " << odd_data.total_mat_energy << "\n";
  std::cout << "     total energy = " << odd_data.total_energy << "\n";
  std::cout << "     relative conservation = "
            << std::abs(odd_data.total_energy - total_energy0) / odd_data.total_energy << "\n";
}

} // namespace odd_driver

//------------------------------------------------------------------------------------------------//
// end of Odd_Functions.cc
//------------------------------------------------------------------------------------------------//
