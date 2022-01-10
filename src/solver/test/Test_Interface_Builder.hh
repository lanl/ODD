//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/test/Test_Interface_Builder.hh
 * \author Mathew Cleveland
 * \brief  Build out simple interface data class for testing. 
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef solver_test_Test_Interface_Builder_hh
#define solver_test_Test_Interface_Builder_hh

#include "solver/Interface_Data.hh"
#include "ds++/dbc.hh"

namespace odd_solver_test {
//================================================================================================//
/*!
 * \brief
 *
 * Build homogeneous single material properties
 *
 */
//================================================================================================//

void Test_Single_Mat_Builder(odd_solver::Interface_Data &iface) {
  Insist(iface.mesh_data.number_of_local_cells > 0, "Must have cells defined on the mesh");
  iface.mat_data.ipcress_filename = "two-mats.ipcress";
  const size_t ncells = iface.mesh_data.number_of_local_cells;
  iface.mat_data.number_of_mats = 1;
  iface.mat_data.problem_matids = {10001};
  iface.mat_data.number_of_cell_mats = std::vector<size_t>(ncells, 1);
  iface.mat_data.cell_mats = std::vector<std::vector<size_t>>(ncells, {0});
  iface.mat_data.cell_mat_vol_frac = std::vector<std::vector<double>>(ncells, {1.0});
  iface.mat_data.cell_mat_temperature = std::vector<std::vector<double>>(ncells, {3.0});
  iface.mat_data.cell_mat_density = std::vector<std::vector<double>>(ncells, {3.0});
  iface.mat_data.cell_mat_specific_heat = std::vector<std::vector<double>>(ncells, {3.0});
  iface.mat_data.cell_rad_eden = std::vector<double>(ncells, 3.0);
  iface.mat_data.cell_velocity = std::vector<std::array<double, 3>>(ncells, {0.0, 0.0, 0.0});
}

//================================================================================================//
/*!
 * \brief
 *
 * Build homogeneous 2 material properties
 *
 */
//================================================================================================//

void Test_Multi_Mat_Builder(odd_solver::Interface_Data &iface) {
  Insist(iface.mesh_data.number_of_local_cells > 0, "Must have cells defined on the mesh");
  iface.mat_data.ipcress_filename = "two-mats.ipcress";
  const size_t ncells = iface.mesh_data.number_of_local_cells;
  iface.mat_data.number_of_mats = 2;
  iface.mat_data.problem_matids = {10001, 10002};
  iface.mat_data.number_of_cell_mats = std::vector<size_t>(ncells, 2);
  iface.mat_data.cell_mats = std::vector<std::vector<size_t>>(ncells, {0, 1});
  iface.mat_data.cell_mat_vol_frac = std::vector<std::vector<double>>(ncells, {0.5, 0.5});
  iface.mat_data.cell_mat_temperature = std::vector<std::vector<double>>(ncells, {3.0, 5.0});
  iface.mat_data.cell_mat_density = std::vector<std::vector<double>>(ncells, {3.0, 5.0});
  iface.mat_data.cell_mat_specific_heat = std::vector<std::vector<double>>(ncells, {3.0, 5.0});
  iface.mat_data.cell_rad_eden = std::vector<double>(ncells, 4.0);
  iface.mat_data.cell_velocity = std::vector<std::array<double, 3>>(ncells, {0.0, 0.0, 0.0});
}

//================================================================================================//
/*!
 * \brief
 *
 * Build out a simple 1D solver interface for testing
 *
 */
//================================================================================================//

void Test_1D_Interface_Builder(odd_solver::Interface_Data &iface) {
  // Define mesh data
  // 2 zones | 0 || 1 | with dx=0.5 dy=0 and dz=0
  iface.mesh_data.domain_decomposed = false;
  iface.mesh_data.number_of_local_cells = 2;
  iface.mesh_data.number_of_global_cells = 2;
  iface.mesh_data.n_dims = 1;
  iface.mesh_data.coord_sys = odd_solver::COORDINATE_SYSTEM::CARTESIAN;
  iface.mesh_data.cell_position = {0.25, 0.0, 0.0, 0.75, 0.0, 0.0};
  iface.mesh_data.cell_size = {0.5, 0.0, 0.0, 0.5, 0.0, 0.0};
  iface.mesh_data.cell_global_id = {0, 1};
  iface.mesh_data.face_types = {
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE};
  iface.mesh_data.next_cell_id = {2, 1, 0, 2};
}

//================================================================================================//
/*!
 * \brief
 *
 * Build out a simple 2D solver interface for testing
 *
 */
//================================================================================================//

void Test_2D_Interface_Builder(odd_solver::Interface_Data &iface) {
  // Define mesh data
  // 4 zones with dx=0.5 dy=0.5 and dz=0.5
  //          ___  ___
  //         | 2 || 3 |
  //          ===  ===
  //         | 0 || 1 |
  //          ---  ---
  iface.mesh_data.domain_decomposed = false;
  iface.mesh_data.number_of_local_cells = 4;
  iface.mesh_data.number_of_global_cells = 4;
  iface.mesh_data.n_dims = 2;
  iface.mesh_data.coord_sys = odd_solver::COORDINATE_SYSTEM::CARTESIAN;
  iface.mesh_data.cell_position = {0.25, 0.25, 0.0, 0.75, 0.25, 0.0,
                                   0.25, 0.75, 0.0, 0.75, 0.75, 0.0};
  iface.mesh_data.cell_size = {0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0, 0.5, 0.5, 0.0};
  iface.mesh_data.cell_global_id = {0, 1, 2, 3};
  iface.mesh_data.face_types = {
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE};
  iface.mesh_data.next_cell_id = {4, 1, 4, 2, 0, 4, 4, 3, 4, 3, 0, 4, 2, 4, 1, 4};
}

void Test_3D_Interface_Builder(odd_solver::Interface_Data &iface) {
  // Define mesh data
  // 8 zones with dx=0.5 dy=0.5 and dz=0.5
  //           front        back
  //          ___  ___    ___  ___
  //         | 2 || 3 |  | 6 || 7 |
  //          ===  ===    ===  ===
  //         | 0 || 1 |  | 4 || 5 |
  //          ---  ---    ---  ---
  iface.mesh_data.domain_decomposed = false;
  iface.mesh_data.number_of_local_cells = 8;
  iface.mesh_data.number_of_global_cells = 8;
  iface.mesh_data.n_dims = 3;
  iface.mesh_data.coord_sys = odd_solver::COORDINATE_SYSTEM::CARTESIAN;
  iface.mesh_data.cell_position = {0.25, 0.25, 0.25, 0.75, 0.25, 0.25, 0.25, 0.75,
                                   0.25, 0.75, 0.75, 0.25, 0.25, 0.25, 0.75, 0.75,
                                   0.25, 0.75, 0.25, 0.75, 0.75, 0.75, 0.75, 0.75};
  iface.mesh_data.cell_size = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5,
                               0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
  iface.mesh_data.cell_global_id = {0, 1, 2, 3, 4, 5, 6, 7};
  iface.mesh_data.face_types = {
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE, // cell 0
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE, // cell 1
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE, // cell 2
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE, // cell 3
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE, // cell 4
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE, // cell 5
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE, // cell 6
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE, // cell 7
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE};
  iface.mesh_data.next_cell_id = {8, 1, 8, 2, 8, 4,  // cell 0
                                  0, 8, 8, 3, 8, 5,  // cell 1
                                  8, 3, 0, 8, 8, 6,  // cell 2
                                  2, 8, 1, 8, 8, 7,  // cell 3
                                  8, 5, 8, 6, 0, 8,  // cell 4
                                  4, 8, 8, 7, 1, 8,  // cell 5
                                  8, 7, 4, 8, 2, 8,  // cell 6
                                  6, 8, 5, 8, 3, 8}; // cell 7
}

} // namespace odd_solver_test

#endif // solver_test_Test_Interface_Builder

//------------------------------------------------------------------------------------------------//
// end of <pkg>/<class>.hh
//------------------------------------------------------------------------------------------------//
