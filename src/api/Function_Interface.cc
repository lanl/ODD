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
#include "solver/Grey_P1_Matrix.hh"
#include "solver/Interface_Data.hh"
#include "solver/MG_Matrix.hh"
#include "solver/MG_P1_Matrix.hh"
#include "c4/global.hh"
#include "ds++/dbc.hh"
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
                                  arg.zonal_data.number_of_cell_mats,
                                  2 * arg.zonal_data.dimensions);

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

  if (arg.control_data.diffusion_method == 0) {
    if (arg.control_data.multigroup == 0) {
      // Initialize a Matrix
      odd_solver::Grey_P1_Matrix matrix(control_data);
      // Initialize Solver Data
      matrix.initialize_solver_data(mesh, mat_data, dt);

      // Build the Matrix Coefficients
      matrix.build_matrix(mesh, dt);

      // Call the solver to calculate the implicit radiation energy density
      matrix.gs_solver(arg.control_data.min_tol, arg.control_data.max_iter,
                       arg.control_data.print == 1);

      // Calculate output data for the implicit radiation energy density vector
      matrix.calculate_output_data(mesh, iface.mat_data, dt, iface.output_data);

    } else {
      // Initialize a Matrix
      odd_solver::MG_P1_Matrix matrix(control_data);

      // Initialize Solver Data
      matrix.initialize_solver_data(mesh, mat_data, dt);

      // Build the Matrix Coefficients
      matrix.build_matrix(mesh, dt);

      // Call the solver to calculate the implicit radiation energy density
      matrix.gs_solver(arg.control_data.min_tol, arg.control_data.max_iter,
                       arg.control_data.print == 1);

      // Calculate output data for the implicit radiation energy density vector
      matrix.calculate_output_data(mesh, iface.mat_data, dt, iface.output_data);
    }
  } else if (arg.control_data.diffusion_method == 1) {
    const bool use_flux_limiter = true;
    if (arg.control_data.multigroup == 0) {
      // Initialize a Matrix
      odd_solver::Grey_Matrix matrix(control_data, use_flux_limiter);

      // Initialize Solver Data
      matrix.initialize_solver_data(mesh, mat_data, dt);

      // Build the Matrix Coefficients
      matrix.build_matrix(mesh, dt);

      // Call the solver to calculate the implicit radiation energy density
      matrix.gs_solver(arg.control_data.min_tol, arg.control_data.max_iter,
                       arg.control_data.print == 1);

      // Calculate output data for the implicit radiation energy density vector
      matrix.calculate_output_data(mesh, iface.mat_data, dt, iface.output_data);

    } else {
      // Initialize a Matrix
      odd_solver::MG_Matrix matrix(control_data, use_flux_limiter);

      // Initialize Solver Data
      matrix.initialize_solver_data(mesh, mat_data, dt);

      // Build the Matrix Coefficients
      matrix.build_matrix(mesh, dt);

      // Call the solver to calculate the implicit radiation energy density
      matrix.gs_solver(arg.control_data.min_tol, arg.control_data.max_iter,
                       arg.control_data.print == 1);

      // Calculate output data for the implicit radiation energy density vector
      matrix.calculate_output_data(mesh, iface.mat_data, dt, iface.output_data);
    }
  } else if (arg.control_data.diffusion_method == 2) {
    // Initialize a Matrix
    const bool use_flux_limiter = false;
    if (arg.control_data.multigroup == 0) {
      // Initialize a Matrix
      odd_solver::Grey_Matrix matrix(control_data, use_flux_limiter);

      // Initialize Solver Data
      matrix.initialize_solver_data(mesh, mat_data, dt);

      // Build the Matrix Coefficients
      matrix.build_matrix(mesh, dt);

      // Call the solver to calculate the implicit radiation energy density
      matrix.gs_solver(arg.control_data.min_tol, arg.control_data.max_iter,
                       arg.control_data.print == 1);

      // Calculate output data for the implicit radiation energy density vector
      matrix.calculate_output_data(mesh, iface.mat_data, dt, iface.output_data);

    } else {
      // Initialize a Matrix
      odd_solver::MG_Matrix matrix(control_data, use_flux_limiter);

      // Initialize Solver Data
      matrix.initialize_solver_data(mesh, mat_data, dt);

      // Build the Matrix Coefficients
      matrix.build_matrix(mesh, dt);

      // Call the solver to calculate the implicit radiation energy density
      matrix.gs_solver(arg.control_data.min_tol, arg.control_data.max_iter,
                       arg.control_data.print == 1);

      // Calculate output data for the implicit radiation energy density vector
      matrix.calculate_output_data(mesh, iface.mat_data, dt, iface.output_data);
    }
  } else {
    Insist(arg.control_data.diffusion_method < 3, "Must specify a valid diffusion method");
  }

  // fill the output arguments
  size_t m_index = 0;
  size_t f_index = 0;
  size_t ngroups = arg.control_data.ngroups;
  for (size_t i = 0; i < mesh_data.number_of_local_cells; i++) {
    arg.output_data.cell_erad[i] = output_data.cell_rad_eden[i];
    for (size_t g = 0; g < ngroups; g++)
      arg.output_data.cell_mg_erad[i * ngroups + g] = output_data.cell_rad_mg_eden[i][g];
    arg.output_data.cell_Trad[i] =
        std::pow(output_data.cell_rad_eden[i] / odd_solver::constants::a, 0.25);
    for (size_t m = 0; m < mat_data.number_of_cell_mats[i]; m++, m_index++)
      arg.output_data.cell_mat_delta_e[m_index] = output_data.cell_mat_dedv[i][m];
    for (size_t f = 0; f < 2 * arg.zonal_data.dimensions; f++, f_index++) {
      arg.output_data.face_flux[f_index] = output_data.face_flux[i][f];
    for (size_t g = 0; g < ngroups; g++)
      arg.output_data.face_mg_flux[f_index * ngroups + g] = output_data.face_mg_flux[i][f][g];
    }
  }
}

void MPI_Initialize(int argc, char *argv[]) {
  // Initialize MPI
  rtt_c4::initialize(argc, argv);
}

//------------------------------------------------------------------------------------------------//
// end of api/Function_Interface.cc
//------------------------------------------------------------------------------------------------//
