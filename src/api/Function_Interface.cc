//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   api/Function_Interface.cc
 * \author Mathew Cleveland
 * \date   November 4th 2021
 * \brief  Contains the C interface function implementations for the ODD solver
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "api/Function_Interface.hh"
#include "solver/Interface_Data.hh"
#include "solver/Opacity_Reader.hh"

//================================================================================================//
/*!
 * \brief Call diffusion solver
 *
 * \param[in] arg interface arguments data
 *
 */
//================================================================================================//
void Odd_Diffusion_Solve(Arguments &arg) {
  // check the argument class to be valid
  arg.control_data.check_arguments();
  arg.zonal_data.check_arguments();
  arg.output_data.check_arguments();

  // Build flat interface (maybe make a builder class for this)
  odd_solver::Interface_Data iface;
  iface.mesh_data.number_of_local_cells = arg.zonal_data.number_of_local_cells;
  iface.mat_data.number_of_mats = arg.zonal_data.number_of_mats;
  iface.mat_data.problem_matids.resize(iface.mat_data.number_of_mats);
  auto matids_ptr = arg.zonal_data.problem_matids;
  for (size_t m = 0; m < iface.mat_data.number_of_mats; m++, matids_ptr++)
    iface.mat_data.problem_matids[m] = *matids_ptr;
  iface.mat_data.number_of_cell_mats.resize(iface.mesh_data.number_of_local_cells);
  auto n_cell_mats_ptr = arg.zonal_data.number_of_cell_mats;
  for (size_t i = 0; i < iface.mesh_data.number_of_local_cells; i++, n_cell_mats_ptr++)
    iface.mat_data.number_of_cell_mats[i] = *n_cell_mats_ptr;
  // unwind the cell-material arrays
  auto cell_mats_ptr = arg.zonal_data.cell_mats;
  auto cell_mat_vol_frac_ptr = arg.zonal_data.cell_mat_vol_frac;
  auto cell_mat_temperature_ptr = arg.zonal_data.cell_mat_temperature;
  auto cell_mat_density_ptr = arg.zonal_data.cell_mat_density;
  iface.mat_data.cell_mats.resize(iface.mesh_data.number_of_local_cells);
  iface.mat_data.cell_mat_vol_frac.resize(iface.mesh_data.number_of_local_cells);
  iface.mat_data.cell_mat_temperature.resize(iface.mesh_data.number_of_local_cells);
  iface.mat_data.cell_mat_density.resize(iface.mesh_data.number_of_local_cells);
  for (size_t i = 0; i < iface.mesh_data.number_of_local_cells; i++) {
    iface.mat_data.cell_mats[i].resize(iface.mat_data.number_of_cell_mats[i]);
    iface.mat_data.cell_mat_vol_frac[i].resize(iface.mat_data.number_of_cell_mats[i]);
    iface.mat_data.cell_mat_temperature[i].resize(iface.mat_data.number_of_cell_mats[i]);
    iface.mat_data.cell_mat_density[i].resize(iface.mat_data.number_of_cell_mats[i]);
    for (size_t m = 0; m < iface.mat_data.number_of_cell_mats[i]; m++) {
      iface.mat_data.cell_mats[i][m] = *cell_mats_ptr;
      iface.mat_data.cell_mat_vol_frac[i][m] = *cell_mat_vol_frac_ptr;
      iface.mat_data.cell_mat_temperature[i][m] = *cell_mat_temperature_ptr;
      iface.mat_data.cell_mat_density[i][m] = *cell_mat_density_ptr;
      cell_mats_ptr++;
      cell_mat_vol_frac_ptr++;
      cell_mat_temperature_ptr++;
      cell_mat_density_ptr++;
    }
  }

  // Open the ipcress file average and set outgoing number of groups
  odd_solver::Opacity_Reader opacity_file(arg.control_data.opacity_file,
                                          iface.mat_data.problem_matids);

  // fill in the opacity data based on volume averages
  auto cell_opacity_ptr = arg.output_data.ave_opacity_data;
  for (size_t i = 0; i < iface.mesh_data.number_of_local_cells; i++) {
    *cell_opacity_ptr = 0.0;
    for (size_t m = 0; m < iface.mat_data.number_of_cell_mats[i]; m++) {
      *cell_opacity_ptr +=
          opacity_file.mat_rosseland_abs_models[m]->getOpacity(
              iface.mat_data.cell_mat_temperature[i][m], iface.mat_data.cell_mat_density[i][m]) *
          iface.mat_data.cell_mat_vol_frac[i][m];
    }
    cell_opacity_ptr++;
  }
}

//------------------------------------------------------------------------------------------------//
// end of api/Function_Interface.cc
//------------------------------------------------------------------------------------------------//
