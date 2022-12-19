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
#include "c4/global.hh"
#include "cdi/CDI.hh"
#include "cdi_analytic/Analytic_EoS.hh"
#include "cdi_analytic/Analytic_Models.hh"
#include "ds++/dbc.hh"
#include <iostream>
#include <string>

namespace odd_driver {

void build_arguments_from_cmd(const std::vector<std::string> &argv, Arguments &args,
                              Odd_Driver_Data &odd_data) {

  // initialize DD mode
  args.zonal_data.domain_decomposed = 0;
  for (size_t i = 0; i < argv.size(); i++) {
    std::string str = argv[i];
    // Parse control data
    if (str == "-if" || str == "-ipcress_file") {
      odd_data.opacity_file = argv[i + 1];
      args.control_data.opacity_file = &odd_data.opacity_file[0];
    }
    if (str == "-s" || str == "-solver") {
      std::string solver = argv[i + 1];
      std::cout << "Solver = " << solver << std::endl;
      if (solver == "P1") {
        args.control_data.diffusion_method = 0;
      } else if (solver == "FLD") {
        args.control_data.diffusion_method = 1;
      } else if (solver == "DIFFUSION") {
        args.control_data.diffusion_method = 2;
      } else {
        Insist(false, "Diffusion option (" + solver + ") not found");
      }
    }
    if (str == "-pf" || str == "-print_frequency")
      odd_data.print_frequency = std::stoi(argv[i + 1]);
    if (str == "-nc" || str == "-n_cycles")
      odd_data.n_cycles = std::stoi(argv[i + 1]);
    if (str == "-dt")
      args.control_data.dt = std::stod(argv[i + 1]);
    if (str == "-dtr")
      odd_data.dt_ramp = std::stod(argv[i + 1]);
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
    // background Mat data
    if (str == "-mid" || str == "-matid")
      odd_data.matid = std::stoi(argv[i + 1]);
    if (str == "-tev" || str == "-temperature")
      odd_data.temperature = std::stod(argv[i + 1]);
    if (str == "-rev" || str == "-rad_temperature")
      odd_data.rad_temperature = std::stod(argv[i + 1]);
    if (str == "-rho" || str == "-density")
      odd_data.density = std::stod(argv[i + 1]);
    if (str == "-aeos" || str == "-analytic_eos") {
      odd_data.specific_heat = std::stod(argv[i + 1]);
      odd_data.specific_heat_Tref = std::stod(argv[i + 2]);
      odd_data.specific_heat_Tpow = std::stod(argv[i + 3]);
      // Convert from jerks/g -> kJ/g
      odd_data.eos =
          std::make_unique<rtt_cdi_analytic::Analytic_EoS>(rtt_cdi_analytic::Analytic_EoS(
              std::make_unique<rtt_cdi_analytic::Polynomial_Specific_Heat_Analytic_EoS_Model>(
                  rtt_cdi_analytic::Polynomial_Specific_Heat_Analytic_EoS_Model(
                      odd_data.specific_heat * 1.0e+6, odd_data.specific_heat_Tref,
                      odd_data.specific_heat_Tpow, 0.0, 0.0, 0.0))));
    }
    if (str == "-vsrc" || str == "-volume_source") {
      odd_data.vol_source_strength = std::stod(argv[i + 1]);
      odd_data.vol_source_duration[0] = std::stod(argv[i + 2]);
      odd_data.vol_source_duration[1] = std::stod(argv[i + 3]);
      odd_data.vol_source_eir_split[0] = std::stod(argv[i + 4]);
      odd_data.vol_source_eir_split[1] = std::stod(argv[i + 5]);
      odd_data.vol_source_eir_split[2] = std::stod(argv[i + 6]);
      odd_data.vol_source_lower_bound[0] = std::stod(argv[i + 7]);
      odd_data.vol_source_lower_bound[1] = std::stod(argv[i + 8]);
      odd_data.vol_source_lower_bound[2] = std::stod(argv[i + 9]);
      odd_data.vol_source_upper_bound[0] = std::stod(argv[i + 10]);
      odd_data.vol_source_upper_bound[1] = std::stod(argv[i + 11]);
      odd_data.vol_source_upper_bound[2] = std::stod(argv[i + 12]);
    }
    if (str == "-bt" || str == "-boundary_temp") {
      odd_data.bnd_temp[0] = std::stod(argv[i + 1]);
      odd_data.bnd_temp[1] = std::stod(argv[i + 2]);
      odd_data.bnd_temp[2] = std::stod(argv[i + 3]);
      odd_data.bnd_temp[3] = std::stod(argv[i + 4]);
      odd_data.bnd_temp[4] = std::stod(argv[i + 5]);
      odd_data.bnd_temp[5] = std::stod(argv[i + 6]);
    }
    if (str == "-rb" || str == "-reflect_bnd") {
      odd_data.reflect_bnd[0] = std::stoi(argv[i + 1]);
      odd_data.reflect_bnd[1] = std::stoi(argv[i + 2]);
      odd_data.reflect_bnd[2] = std::stoi(argv[i + 3]);
      odd_data.reflect_bnd[3] = std::stoi(argv[i + 4]);
      odd_data.reflect_bnd[4] = std::stoi(argv[i + 5]);
      odd_data.reflect_bnd[5] = std::stoi(argv[i + 6]);
    }
    if (str == "-dd" || str == "-domain_decomposed") {
      args.zonal_data.domain_decomposed = 1;
      odd_data.domain_decomposed = true;
    }
    // Specify problem regions
    if (str == "-br" || str == "-block_region") {
      odd_driver::Odd_Region block_region;
      block_region.block = true;
      block_region.matid = std::stoi(argv[i + 1]);
      block_region.temperature = std::stod(argv[i + 2]);
      block_region.rad_temperature = std::stod(argv[i + 3]);
      block_region.density = std::stod(argv[i + 4]);
      block_region.specific_heat = std::stod(argv[i + 5]);
      block_region.specific_heat_Tref = std::stod(argv[i + 6]);
      block_region.specific_heat_Tpow = std::stod(argv[i + 7]);
      // Convert from jerks/g -> kJ/g
      block_region.eos =
          std::make_unique<rtt_cdi_analytic::Analytic_EoS>(rtt_cdi_analytic::Analytic_EoS(
              std::make_unique<rtt_cdi_analytic::Polynomial_Specific_Heat_Analytic_EoS_Model>(
                  rtt_cdi_analytic::Polynomial_Specific_Heat_Analytic_EoS_Model(
                      block_region.specific_heat * 1.0e+6, block_region.specific_heat_Tref,
                      block_region.specific_heat_Tpow, 0.0, 0.0, 0.0))));
      block_region.block_p0[0] = std::stod(argv[i + 8]);
      block_region.block_p0[1] = std::stod(argv[i + 9]);
      block_region.block_p0[2] = std::stod(argv[i + 10]);
      block_region.block_p1[0] = std::stod(argv[i + 11]);
      block_region.block_p1[1] = std::stod(argv[i + 12]);
      block_region.block_p1[2] = std::stod(argv[i + 13]);
      block_region.region_number = odd_data.regions.size() + 1;
      odd_data.regions.push_back(std::move(block_region));
    }
    if (str == "-sr" || str == "-sphere_region") {
      Odd_Region sphere_region;
      sphere_region.sphere = true;
      sphere_region.matid = std::stoi(argv[i + 1]);
      sphere_region.temperature = std::stod(argv[i + 2]);
      sphere_region.rad_temperature = std::stod(argv[i + 3]);
      sphere_region.density = std::stod(argv[i + 4]);
      sphere_region.specific_heat = std::stod(argv[i + 5]);
      sphere_region.specific_heat_Tref = std::stod(argv[i + 6]);
      sphere_region.specific_heat_Tpow = std::stod(argv[i + 7]);
      // Convert from jerks/g -> kJ/g
      sphere_region.eos =
          std::make_unique<rtt_cdi_analytic::Analytic_EoS>(rtt_cdi_analytic::Analytic_EoS(
              std::make_unique<rtt_cdi_analytic::Polynomial_Specific_Heat_Analytic_EoS_Model>(
                  rtt_cdi_analytic::Polynomial_Specific_Heat_Analytic_EoS_Model(
                      sphere_region.specific_heat * 1.0e+6, sphere_region.specific_heat_Tref,
                      sphere_region.specific_heat_Tpow, 0.0, 0.0, 0.0))));
      sphere_region.sphere_radius = std::stod(argv[i + 8]);
      sphere_region.sphere_center[0] = std::stod(argv[i + 9]);
      sphere_region.sphere_center[1] = std::stod(argv[i + 10]);
      sphere_region.sphere_center[2] = std::stod(argv[i + 11]);
      sphere_region.region_number = odd_data.regions.size() + 1;
      odd_data.regions.push_back(std::move(sphere_region));
    }
  }
  // Assign remaining control data
  args.control_data.bnd_temp = &odd_data.bnd_temp[0];
  args.control_data.reflect_bnd = &odd_data.reflect_bnd[0];

  // Build mesh arguments
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
                    ? odd_data.mesh_n_cells[1]
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

  // simple decomposition
  if (args.zonal_data.domain_decomposed == 1) {
    const size_t equal_ncells = std::max((ncells / static_cast<size_t>(rtt_c4::nodes())), 1UL);
    size_t global_ncells = equal_ncells;
    rtt_c4::global_sum(global_ncells);
    Insist(!(global_ncells > ncells),
           "The number of parallel global cells can not be more then the number of serial cells");

    size_t local_ncells = equal_ncells;
    if (rtt_c4::node() == (rtt_c4::nodes() - 1)) {
      // add the remaining ranks to the last proc
      size_t remainder = ncells - static_cast<size_t>(rtt_c4::nodes()) * equal_ncells;
      local_ncells += remainder;
    }
    global_ncells = local_ncells;
    rtt_c4::global_sum(global_ncells);
    Insist(global_ncells == ncells, "Global cells does not match mesh specification");

    const size_t cell_id_initial = static_cast<size_t>(rtt_c4::node()) * equal_ncells;
    const size_t cell_id_final = cell_id_initial + local_ncells;
    args.zonal_data.number_of_local_cells = local_ncells;

    // make a local copy of my partition data
    std::vector<double> local_cell_position(odd_data.cell_position.begin() + cell_id_initial * 3,
                                            odd_data.cell_position.begin() + cell_id_final * 3);
    std::vector<double> local_cell_size(odd_data.cell_size.begin() + cell_id_initial * 3,
                                        odd_data.cell_size.begin() + cell_id_final * 3);
    std::vector<size_t> local_cell_id(odd_data.cell_global_id.begin() + cell_id_initial,
                                      odd_data.cell_global_id.begin() + cell_id_final);
    std::vector<size_t> local_face_type(
        odd_data.face_type.begin() + cell_id_initial * nfaces_per_cell,
        odd_data.face_type.begin() + cell_id_final * nfaces_per_cell);
    std::vector<size_t> local_next_cell_id(
        odd_data.next_cell_id.begin() + cell_id_initial * nfaces_per_cell,
        odd_data.next_cell_id.begin() + cell_id_final * nfaces_per_cell);

    // Initialize ghost map from local global data
    // ghost_map[global_id] = {proc, ghost_index}
    std::map<size_t, std::array<size_t, 2>> local_ghost_map;
    for (size_t i = 0; i < local_face_type.size(); i++) {
      const size_t gid = local_next_cell_id[i];
      if ((gid < cell_id_initial || gid >= cell_id_final) && gid < ncells) {
        // setup the local ghost map
        const size_t proc = std::min((gid / equal_ncells), static_cast<size_t>(rtt_c4::nodes()));
        local_ghost_map[gid] = {proc, 0UL};
      }
    }
    // now set the ghost data;
    std::vector<size_t> local_ghost_ids(local_ghost_map.size(), 0UL);
    std::vector<size_t> local_ghost_procs(local_ghost_map.size(), 0UL);
    size_t local_gid = 0;
    for (auto &map : local_ghost_map) {
      local_ghost_ids[local_gid] = map.first;
      local_ghost_procs[local_gid] = map.second[0];
      map.second[1] = local_gid;
      local_gid++;
    }

    // Alright now we can clean up the local face types and next cell ids
    for (size_t i = 0; i < local_face_type.size(); i++) {
      const size_t gid = local_next_cell_id[i];
      local_next_cell_id[i] = gid - cell_id_initial;
      if (gid < cell_id_initial || gid >= cell_id_final) {
        // reset the face type (GHOST OR BOUNDRAY=ncells)
        if (gid < ncells) {
          Check(local_face_type[i] == 0);
          // set as ghost
          local_face_type[i] = 2;
          // give it the ghost index;
          local_next_cell_id[i] = local_ghost_map[gid][1];
        } else {
          // fix local
          Check(local_face_type[i] == 1);
          local_next_cell_id[i] = local_ncells;
        }
      }
    }

    // Now reset the ODD data
    args.zonal_data.number_of_local_cells = local_ncells;
    args.zonal_data.number_of_ghost_cells = local_ghost_ids.size();
    args.zonal_data.number_of_global_cells = global_ncells;
    odd_data.cell_position = local_cell_position;
    odd_data.cell_size = local_cell_size;
    odd_data.cell_global_id = local_cell_id;
    odd_data.face_type = local_face_type;
    odd_data.next_cell_id = local_next_cell_id;
    odd_data.ghost_cell_global_id = local_ghost_ids;
    args.zonal_data.ghost_cell_global_id = &odd_data.ghost_cell_global_id[0];
    odd_data.ghost_cell_proc = local_ghost_procs;
    args.zonal_data.ghost_cell_proc = &odd_data.ghost_cell_proc[0];

    // rest ncells for the material population
    ncells = local_ncells;
  }

  // Allocate and assign the material data
  args.zonal_data.number_of_mats = odd_data.regions.size() + 1;
  odd_data.problem_matids = std::vector<size_t>(args.zonal_data.number_of_mats, odd_data.matid);
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
  odd_data.cell_velocity = std::vector<double>(ncells * 3, 0.0);
  args.zonal_data.cell_velocity = &odd_data.cell_velocity[0];
  odd_data.cell_erad = std::vector<double>(ncells, odd_solver::constants::a *
                                                       std::pow(odd_data.rad_temperature, 4.0));
  args.zonal_data.cell_erad = &odd_data.cell_erad[0];
  odd_data.face_flux = std::vector<double>(ncells * nfaces_per_cell, 0.0);
  args.zonal_data.face_flux = &odd_data.face_flux[0];

  // Allocate mat source arrays
  odd_data.cell_mat_electron_source = std::vector<double>(ncells, 0.0);
  args.zonal_data.cell_mat_electron_source = &odd_data.cell_mat_electron_source[0];
  odd_data.cell_rad_source = std::vector<double>(ncells, 0.0);
  args.zonal_data.cell_rad_source = &odd_data.cell_rad_source[0];

  // Output data
  odd_data.output_cell_erad = std::vector<double>(
      ncells, odd_solver::constants::a * std::pow(odd_data.rad_temperature, 4.0));
  args.output_data.cell_erad = &odd_data.output_cell_erad[0];
  odd_data.output_cell_Trad = std::vector<double>(ncells, odd_data.rad_temperature);
  args.output_data.cell_Trad = &odd_data.output_cell_Trad[0];
  odd_data.output_cell_mat_delta_e = std::vector<double>(ncells, 0.0);
  args.output_data.cell_mat_delta_e = &odd_data.output_cell_mat_delta_e[0];
  odd_data.output_face_flux = std::vector<double>(ncells * nfaces_per_cell, 0.0);
  args.output_data.face_flux = &odd_data.output_face_flux[0];

  // Query the EOS in kJ/g
  odd_data.cell_mat_specific_heat = odd_data.eos->getElectronHeatCapacity(
      odd_data.cell_mat_temperature, odd_data.cell_mat_density);
  args.zonal_data.cell_mat_specific_heat = &odd_data.cell_mat_specific_heat[0];
  odd_data.cell_mat_energy_density = odd_data.eos->getSpecificElectronInternalEnergy(
      odd_data.cell_mat_temperature, odd_data.cell_mat_density);

  for (size_t r = 0; r < odd_data.regions.size(); r++)
    odd_data.problem_matids[r + 1] = odd_data.regions[r].matid;
  size_t mat_index = 0;
  odd_data.total_rad_energy = 0.0;
  odd_data.total_mat_energy = 0.0;
  for (size_t i = 0; i < args.zonal_data.number_of_local_cells; i++) {
    // Paint on regions in the order they are parsed in the input
    // (still single material). Calculate
    size_t eos_index = 0;
    for (const auto &region : odd_data.regions) {
      bool valid_cell =
          region.block
              ? odd_data.cell_position[i * 3 + 0] <= region.block_p1[0] &&
                    odd_data.cell_position[i * 3 + 0] >= region.block_p0[0] &&
                    odd_data.cell_position[i * 3 + 1] <= region.block_p1[1] &&
                    odd_data.cell_position[i * 3 + 1] >= region.block_p0[1] &&
                    odd_data.cell_position[i * 3 + 2] <= region.block_p1[2] &&
                    odd_data.cell_position[i * 3 + 2] >= region.block_p0[2]
              : std::sqrt(
                    std::pow(odd_data.cell_position[i * 3 + 0] - region.sphere_center[0], 2.0) +
                    std::pow(odd_data.cell_position[i * 3 + 1] - region.sphere_center[1], 2.0) +
                    std::pow(odd_data.cell_position[i * 3 + 2] - region.sphere_center[2], 2.0)) <=
                    region.sphere_radius;
      // This indexing will not work for multimaterial definitions
      if (valid_cell) {
        Insist(args.zonal_data.number_of_cell_mats[i] == 1,
               "Block and Sphere Regions do not support multimaterial yet");
        odd_data.cell_mats[i] = region.region_number;
        odd_data.cell_mat_temperature[i] = region.temperature;
        odd_data.cell_mat_density[i] = region.density;
        odd_data.cell_mat_specific_heat[i] =
            region.eos->getElectronHeatCapacity(region.temperature, region.density);
        odd_data.cell_mat_energy_density[i] =
            region.eos->getSpecificElectronInternalEnergy(region.temperature, region.density);
        odd_data.cell_erad[i] = odd_solver::constants::a * std::pow(region.rad_temperature, 4.0);
        odd_data.output_cell_erad[i] =
            odd_solver::constants::a * std::pow(region.rad_temperature, 4.0);
        odd_data.output_cell_Trad[i] = region.rad_temperature;
      }
      eos_index++;
    }
    // Convert from kJ/g -> jerks/g
    double volume = 1.0;
    for (size_t d = 0; d < args.zonal_data.dimensions; d++)
      volume *= odd_data.cell_size[i * 3 + d];
    odd_data.total_rad_energy += odd_data.cell_erad[i] * volume;
    for (size_t m = 0; m < args.zonal_data.number_of_cell_mats[i]; m++, mat_index++) {
      // Convert from kJ/g -> jerks/g
      odd_data.cell_mat_specific_heat[mat_index] *= 1.0e-6;
      odd_data.cell_mat_energy_density[mat_index] *= 1.0e-6 * odd_data.cell_mat_density[mat_index];
      odd_data.total_mat_energy += odd_data.cell_mat_energy_density[mat_index] *
                                   odd_data.cell_mat_vol_frac[mat_index] * volume;
    }
  }
  if (odd_data.domain_decomposed) {
    rtt_c4::global_sum(odd_data.total_mat_energy);
    rtt_c4::global_sum(odd_data.total_rad_energy);
  }
  odd_data.total_energy = odd_data.total_mat_energy + odd_data.total_rad_energy;
}

void energy_update(Arguments &args, Odd_Driver_Data &odd_data, bool print_info) {

  if (print_info && rtt_c4::node() == 0) {
    std::cout << " CONSERVATION DATA \n";
    std::cout << "     total rad_energy0 = " << odd_data.total_rad_energy << "\n";
    std::cout << "     total mat_energy0 = " << odd_data.total_mat_energy << "\n";
    std::cout << "     total energy0 = " << odd_data.total_energy << "\n";
  }
  const double total_energy0 = odd_data.total_energy;
  size_t mat_index = 0;
  size_t face_index = 0;
  const size_t nfaces_per_cell = 2 * args.zonal_data.dimensions;
  odd_data.total_rad_energy = 0.0;
  odd_data.total_mat_energy = 0.0;
  odd_data.total_source_energy = 0.0;
  for (size_t i = 0; i < args.zonal_data.number_of_local_cells; i++) {
    double volume = 1.0;
    for (size_t d = 0; d < args.zonal_data.dimensions; d++)
      volume *= odd_data.cell_size[i * 3 + d];
    odd_data.cell_erad[i] = args.output_data.cell_erad[i];
    odd_data.total_rad_energy += odd_data.cell_erad[i] * volume;
    odd_data.total_source_energy += odd_data.cell_rad_source[i] * volume;
    for (size_t m = 0; m < args.zonal_data.number_of_cell_mats[i]; m++, mat_index++) {
      // fetch the correct eos for the current zone material
      rtt_cdi::EoS &eos = args.zonal_data.cell_mats[m] > 0
                              ? *odd_data.regions[args.zonal_data.cell_mats[m] - 1].eos
                              : *odd_data.eos;

      odd_data.cell_mat_energy_density[mat_index] += args.output_data.cell_mat_delta_e[mat_index];
      odd_data.total_mat_energy += odd_data.cell_mat_energy_density[mat_index] *
                                   odd_data.cell_mat_vol_frac[mat_index] * volume;
      odd_data.total_source_energy += odd_data.cell_mat_electron_source[mat_index] *
                                      odd_data.cell_mat_vol_frac[mat_index] * volume;

      // Convert from jerks/g -> kJ/g
      odd_data.cell_mat_temperature[mat_index] =
          eos.getElectronTemperature(odd_data.cell_mat_density[mat_index],
                                     odd_data.cell_mat_energy_density[mat_index] * 1.0e6 /
                                         odd_data.cell_mat_density[mat_index],
                                     odd_data.cell_mat_temperature[mat_index]);
    }
    for (size_t f = 0; f < nfaces_per_cell; f++, face_index++) {
      odd_data.face_flux[face_index] = args.output_data.face_flux[face_index];
    }
  }
  if (odd_data.domain_decomposed) {
    rtt_c4::global_sum(odd_data.total_mat_energy);
    rtt_c4::global_sum(odd_data.total_rad_energy);
    rtt_c4::global_sum(odd_data.total_source_energy);
  }
  odd_data.total_energy = odd_data.total_mat_energy + odd_data.total_rad_energy;

  if (print_info && rtt_c4::node() == 0) {
    std::cout << "     total source_energy = " << odd_data.total_source_energy << "\n";
    std::cout << "     total rad_energy = " << odd_data.total_rad_energy << "\n";
    std::cout << "     total mat_energy = " << odd_data.total_mat_energy << "\n";
    std::cout << "     total energy = " << odd_data.total_energy << "\n";
    std::cout << "     relative conservation = "
              << std::abs(odd_data.total_energy - (total_energy0 + odd_data.total_source_energy)) /
                     odd_data.total_energy
              << "\n";
  }
}

void update_source(Arguments &args, Odd_Driver_Data &odd_data, const double time) {
  const double normal = odd_data.vol_source_eir_split[0] + odd_data.vol_source_eir_split[1] +
                        odd_data.vol_source_eir_split[2];
  const double eratio = normal > 0 ? odd_data.vol_source_eir_split[0] / normal : 0.0;
  const double rratio = normal > 0 ? odd_data.vol_source_eir_split[2] / normal : 0.0;
  Insist(rtt_dsxx::soft_equiv(odd_data.vol_source_eir_split[1], 0.0),
         "Ion sources are not currently supported so eir_split for ions must be zero");
  size_t mat_index = 0;
  for (size_t i = 0; i < args.zonal_data.number_of_local_cells; i++) {
    const size_t n_mats = odd_data.number_of_cell_mats[i];
    const std::array<double, 3> position{odd_data.cell_position[i * 3],
                                         odd_data.cell_position[i * 3 + 1],
                                         odd_data.cell_position[i * 3 + 2]};
    if (position[0] <= odd_data.vol_source_upper_bound[0] &&
        position[0] >= odd_data.vol_source_lower_bound[0] &&
        position[1] <= odd_data.vol_source_upper_bound[1] &&
        position[1] >= odd_data.vol_source_lower_bound[1] &&
        position[2] <= odd_data.vol_source_upper_bound[2] &&
        position[2] >= odd_data.vol_source_lower_bound[2] &&
        time >= odd_data.vol_source_duration[0] && time <= odd_data.vol_source_duration[1]) {
      odd_data.cell_rad_source[i] = rratio * odd_data.vol_source_strength * args.control_data.dt;
      for (size_t mat = 0; mat < n_mats; mat++, mat_index++) {
        odd_data.cell_mat_electron_source[mat_index] =
            eratio * odd_data.vol_source_strength * args.control_data.dt;
      }
    } else {
      odd_data.cell_rad_source[i] = 0.0;
      for (size_t mat = 0; mat < n_mats; mat++, mat_index++) {
        odd_data.cell_mat_electron_source[mat_index] = 0.0;
      }
    }
  }
}

} // namespace odd_driver

//------------------------------------------------------------------------------------------------//
// end of Odd_Functions.cc
//------------------------------------------------------------------------------------------------//
