//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   api/Interface_Builder.cc
 * \author Mathew Cleveland
 * \date   January 13th 2022
 * \brief  Function to manage the Solver:Interface_Data and API data interface
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Interface_Builder.hh"
#include <iostream>

namespace odd_api {

odd_solver::Interface_Data build_interface_data(const Arguments &arg) {
  odd_solver::Interface_Data iface;

  // Populate control data
  for (size_t f = 0; f < 6; f++) {
    iface.control_data.reflect_bnd[f] = (arg.control_data.reflect_bnd[f] == 1);
    iface.control_data.bnd_temp[f] = arg.control_data.bnd_temp[f];
  }
  iface.control_data.correction = (arg.control_data.correction == 1);
  iface.control_data.multigroup = (arg.control_data.multigroup == 1);

  //construct the zonal data
  arg.zonal_data.check_arguments();
  auto &mesh_data = iface.mesh_data;
  mesh_data.domain_decomposed = arg.zonal_data.domain_decomposed;
  mesh_data.n_dims = arg.zonal_data.dimensions;
  mesh_data.coord_sys = arg.zonal_data.coord_sys;
  const size_t ncells = arg.zonal_data.number_of_local_cells;
  mesh_data.number_of_local_cells = ncells;
  const size_t ngcells = arg.zonal_data.number_of_ghost_cells;
  mesh_data.number_of_ghost_cells = ngcells;
  mesh_data.number_of_global_cells = arg.zonal_data.number_of_global_cells;
  mesh_data.cell_position =
      std::vector<double>(arg.zonal_data.cell_position, arg.zonal_data.cell_position + ncells * 3);
  mesh_data.cell_size =
      std::vector<double>(arg.zonal_data.cell_size, arg.zonal_data.cell_size + ncells * 3);
  mesh_data.cell_global_id =
      std::vector<size_t>(arg.zonal_data.cell_global_id, arg.zonal_data.cell_global_id + ncells);
  if (mesh_data.domain_decomposed != 0) {
    mesh_data.ghost_cell_global_id = std::vector<size_t>(
        arg.zonal_data.ghost_cell_global_id, arg.zonal_data.ghost_cell_global_id + ngcells);
    mesh_data.ghost_cell_proc = std::vector<size_t>(arg.zonal_data.ghost_cell_proc,
                                                    arg.zonal_data.ghost_cell_proc + ngcells);
  }
  const size_t nfaces = 2 * mesh_data.n_dims;
  mesh_data.face_types =
      std::vector<size_t>(arg.zonal_data.face_type, arg.zonal_data.face_type + ncells * nfaces);
  mesh_data.next_cell_id = std::vector<size_t>(arg.zonal_data.next_cell_id,
                                               arg.zonal_data.next_cell_id + ncells * nfaces);

  // Populate the material data
  auto &mat_data = iface.mat_data;
  mat_data.ipcress_filename = arg.control_data.opacity_file;
  mat_data.number_of_mats = arg.zonal_data.number_of_mats;
  mat_data.problem_matids = std::vector<size_t>(
      arg.zonal_data.problem_matids, arg.zonal_data.problem_matids + arg.zonal_data.number_of_mats);
  mat_data.number_of_cell_mats = std::vector<size_t>(arg.zonal_data.number_of_cell_mats,
                                                     arg.zonal_data.number_of_cell_mats + ncells);
  // unwind the cell-material arrays
  // resize cell wise arrays
  mat_data.cell_mats.resize(ncells);
  mat_data.cell_mat_vol_frac.resize(iface.mesh_data.number_of_local_cells);
  mat_data.cell_mat_temperature.resize(iface.mesh_data.number_of_local_cells);
  mat_data.cell_mat_density.resize(iface.mesh_data.number_of_local_cells);
  mat_data.cell_mat_specific_heat.resize(iface.mesh_data.number_of_local_cells);
  mat_data.cell_mat_electron_source.resize(iface.mesh_data.number_of_local_cells);
  size_t mat_index = 0;
  for (size_t i = 0; i < iface.mesh_data.number_of_local_cells; i++) {
    // resize cell material arrays
    mat_data.cell_mats[i].resize(mat_data.number_of_cell_mats[i]);
    mat_data.cell_mat_vol_frac[i].resize(mat_data.number_of_cell_mats[i]);
    mat_data.cell_mat_temperature[i].resize(mat_data.number_of_cell_mats[i]);
    mat_data.cell_mat_density[i].resize(mat_data.number_of_cell_mats[i]);
    mat_data.cell_mat_specific_heat[i].resize(mat_data.number_of_cell_mats[i]);
    mat_data.cell_mat_electron_source[i].resize(mat_data.number_of_cell_mats[i]);
    // populate arrays
    for (size_t m = 0; m < mat_data.number_of_cell_mats[i]; m++, mat_index++) {
      mat_data.cell_mats[i][m] = arg.zonal_data.cell_mats[mat_index];
      mat_data.cell_mat_vol_frac[i][m] = arg.zonal_data.cell_mat_vol_frac[mat_index];
      mat_data.cell_mat_temperature[i][m] = arg.zonal_data.cell_mat_temperature[mat_index];
      mat_data.cell_mat_density[i][m] = arg.zonal_data.cell_mat_density[mat_index];
      mat_data.cell_mat_specific_heat[i][m] = arg.zonal_data.cell_mat_specific_heat[mat_index];
      mat_data.cell_mat_electron_source[i][m] = arg.zonal_data.cell_mat_electron_source[mat_index];
    }
  }

  // cell material data
  mat_data.cell_rad_eden =
      std::vector<double>(arg.zonal_data.cell_erad, arg.zonal_data.cell_erad + ncells);
  mat_data.cell_rad_source =
      std::vector<double>(arg.zonal_data.cell_rad_source, arg.zonal_data.cell_rad_source + ncells);
  mat_data.cell_velocity.resize(ncells, {32.0, 32.0, 32.0});
  size_t v_index = 0;
  for (auto &cell_v : mat_data.cell_velocity) {
    for (auto &v : cell_v) {
      v = arg.zonal_data.cell_velocity[v_index];
      v_index++;
    }
  }
  mat_data.face_flux = std::vector<std::vector<double>>(ncells, std::vector<double>(nfaces, 32.0));
  // fill in the face_flux data for p1 solver
  size_t f_index = 0;
  for (auto &flux : mat_data.face_flux) {
    for (auto &face_flux : flux) {
      face_flux = arg.zonal_data.face_flux[f_index];
      f_index++;
    }
  }

  // build output_data
  auto &output_data = iface.output_data;
  output_data.cell_rad_eden = std::vector<double>(ncells, 0.0);
  output_data.cell_mat_dedv = std::vector<std::vector<double>>(ncells);
  for (size_t i = 0; i < iface.mesh_data.number_of_local_cells; i++)
    output_data.cell_mat_dedv[i] = std::vector<double>(mat_data.number_of_cell_mats[i], 0.0);
  output_data.face_flux =
      std::vector<std::vector<double>>(ncells, std::vector<double>(nfaces, 0.0));

  return iface;
}

} // namespace odd_api

//------------------------------------------------------------------------------------------------//
// end of Interface_Data.cc
//------------------------------------------------------------------------------------------------//
