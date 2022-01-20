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
#include "Interface_Builder.hh"
#include "solver/Constants.hh"
#include "solver/Grey_Matrix.hh"
#include "solver/Interface_Data.hh"
#include <cmath>
#include <iostream>

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
  arg.output_data.check_arguments(arg.zonal_data.number_of_local_cells,
                                  arg.zonal_data.number_of_cell_mats);

  // populate the interface data
  odd_solver::Interface_Data iface = odd_api::build_interface_data(arg);
  // Data accessors
  const auto &control_data = iface.control_data;
  const auto &mesh_data = iface.mesh_data;
  const auto &mat_data = iface.mat_data;
  auto &output_data = iface.output_data;
  const auto dt = arg.control_data.dt;

  // Build a Mesh
  odd_solver::Orthogonal_Mesh mesh(mesh_data);

  // Initialize a Matrix
  odd_solver::Grey_Matrix matrix(control_data);

  // Initialize Solver Data
  matrix.initialize_solver_data(mesh, mat_data, dt);

  // Build the Matrix Coefficients
  matrix.build_matrix(mesh, dt);

  // Call the solver to calculate the implicit radiation energy density
  matrix.gs_solver(arg.control_data.min_tol, arg.control_data.max_iter);

  // Calculate output data for the implicit radiation energy density vector
  matrix.calculate_output_data(iface.mat_data, dt, iface.output_data);

  // fill the output arguments
  size_t m_index = 0;
  for (size_t i = 0; i < mesh_data.number_of_local_cells; i++) {
    arg.output_data.cell_erad[i] = output_data.cell_rad_eden[i];
    arg.output_data.cell_Trad[i] =
        std::pow(output_data.cell_rad_eden[i] / odd_solver::constants::a, 0.25);
    for (size_t m = 0; m < mat_data.number_of_cell_mats[i]; m++, m_index++)
      arg.output_data.cell_mat_delta_e[m_index] = output_data.cell_mat_dedv[i][m];
  }
}

//------------------------------------------------------------------------------------------------//
// end of api/Function_Interface.cc
//------------------------------------------------------------------------------------------------//
