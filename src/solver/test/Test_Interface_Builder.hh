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
#include "c4/global.hh"
#include "ds++/dbc.hh"

namespace odd_solver_test {

//================================================================================================//
/*!
 * \brief
 *
 * Initialize output data
 *
 */
//================================================================================================//

void Test_Output_Builder(odd_solver::Interface_Data &iface) {
  Insist(iface.mesh_data.number_of_local_cells > 0, "Must have cells defined on the mesh");
  iface.mat_data.ipcress_filename = "two-mats.ipcress";
  const size_t ncells = iface.mesh_data.number_of_local_cells;
  iface.output_data.cell_rad_eden = std::vector<double>(ncells, 0.0);
  iface.output_data.cell_mat_dedv = std::vector<std::vector<double>>(ncells);
  for (size_t cell = 0; cell < ncells; cell++)
    iface.output_data.cell_mat_dedv[cell].resize(iface.mat_data.number_of_cell_mats[cell], 0.0);
  iface.output_data.face_flux = std::vector<std::vector<double>>(
      ncells, std::vector<double>(iface.mesh_data.n_dims * 2, 0.0));
}

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
  iface.mat_data.face_flux = std::vector<std::vector<double>>(
      ncells, std::vector<double>(iface.mesh_data.n_dims * 2, 0.0));
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

void Test_1D_Interface_Builder(odd_solver::Interface_Data &iface, const bool dd = false) {
  size_t n_global_cells = 2;
  size_t n_local_cells = n_global_cells;
  if (dd) {
    Require(rtt_c4::nodes() < 3);
    n_local_cells = n_global_cells / static_cast<size_t>(rtt_c4::nodes());
  }
  std::vector<double> global_cell_positions = {0.25, 0.0, 0.0, 0.75, 0.0, 0.0};
  std::vector<double> global_cell_size = {0.5, 0.0, 0.0, 0.5, 0.0, 0.0};
  std::vector<size_t> global_cell_id = {0, 1};
  std::vector<size_t> global_next_cell_id = {2, 1, 0, 2};
  // Define mesh data
  // 2 zones | 0 || 1 | with dx=0.5 dy=0 and dz=0
  iface.mesh_data.domain_decomposed = dd;
  iface.mesh_data.number_of_local_cells = n_local_cells;
  iface.mesh_data.number_of_global_cells = n_global_cells;
  iface.mesh_data.n_dims = 1;
  iface.mesh_data.coord_sys = odd_solver::COORDINATE_SYSTEM::CARTESIAN;
  if (!dd) {
    iface.mesh_data.cell_position = global_cell_positions;
    iface.mesh_data.cell_size = global_cell_size;
    iface.mesh_data.cell_global_id = global_cell_id;
    iface.mesh_data.next_cell_id = global_next_cell_id;
    iface.mesh_data.face_types = {
        odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
        odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE};
  } else {
    iface.mesh_data.cell_position.assign(
        global_cell_positions.begin() + rtt_c4::node() * n_local_cells * 3,
        global_cell_positions.begin() + rtt_c4::node() * n_local_cells * 3 + n_local_cells * 3);
    iface.mesh_data.cell_size.assign(global_cell_size.begin() + rtt_c4::node() * n_local_cells * 3,
                                     global_cell_size.begin() + rtt_c4::node() * n_local_cells * 3 +
                                         n_local_cells * 3);
    iface.mesh_data.cell_global_id.assign(global_cell_id.begin() + rtt_c4::node() * n_local_cells,
                                          global_cell_id.begin() + rtt_c4::node() * n_local_cells +
                                              n_local_cells);
    iface.mesh_data.number_of_ghost_cells = 1;
    if (rtt_c4::node() == 0) {
      iface.mesh_data.face_types = {odd_solver::FACE_TYPE::BOUNDARY_FACE,
                                    odd_solver::FACE_TYPE::GHOST_FACE};
      iface.mesh_data.ghost_cell_global_id = {1};
      iface.mesh_data.ghost_cell_proc = {1};
      iface.mesh_data.next_cell_id = {2, 0};
    } else {
      iface.mesh_data.face_types = {odd_solver::FACE_TYPE::GHOST_FACE,
                                    odd_solver::FACE_TYPE::BOUNDARY_FACE};
      iface.mesh_data.ghost_cell_global_id = {0};
      iface.mesh_data.ghost_cell_proc = {0};
      iface.mesh_data.next_cell_id = {0, 2};
    }
  }
}

//================================================================================================//
/*!
 * \brief
 *
 * Build out a simple 2D solver interface for testing
 *
 */
//================================================================================================//

void Test_2D_Interface_Builder(odd_solver::Interface_Data &iface, const bool dd = false) {
  // Define mesh data
  // 4 zones with dx=0.5 dy=0.5 and dz=0.5
  //          ___  ___
  //         | 2 || 3 |
  //          ===  ===
  //         | 0 || 1 |
  //          ---  ---
  // if dd -> cells_node_0={0} cell_node_1={1] cell_node_2={2,3}
  size_t n_global_cells = 4;
  size_t n_local_cells = n_global_cells;
  if (dd) {
    Insist(rtt_c4::nodes() == 3, "2D DD test mesh only supports 3 ranks");
    if (rtt_c4::node() < 2)
      n_local_cells = 1;
    else
      n_local_cells = 2;
  }
  iface.mesh_data.domain_decomposed = dd;
  iface.mesh_data.number_of_local_cells = n_local_cells;
  iface.mesh_data.number_of_global_cells = n_global_cells;
  iface.mesh_data.n_dims = 2;
  iface.mesh_data.coord_sys = odd_solver::COORDINATE_SYSTEM::CARTESIAN;
  if (!dd) {
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
  } else {
    if (rtt_c4::node() == 0) {
      iface.mesh_data.cell_position = {0.25, 0.25, 0.0};
      iface.mesh_data.cell_size = {0.5, 0.5, 0.0};
      iface.mesh_data.cell_global_id = {0};
      iface.mesh_data.face_types = {
          odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::GHOST_FACE,
          odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::GHOST_FACE};
      iface.mesh_data.number_of_ghost_cells = 2;
      iface.mesh_data.next_cell_id = {4, 0, 4, 1};
      iface.mesh_data.ghost_cell_global_id = {1, 2};
      iface.mesh_data.ghost_cell_proc = {1, 2};

    } else if (rtt_c4::node() == 1) {
      iface.mesh_data.cell_position = {0.75, 0.25, 0.0};
      iface.mesh_data.cell_size = {0.5, 0.5, 0.0};
      iface.mesh_data.cell_global_id = {1};
      iface.mesh_data.face_types = {
          odd_solver::FACE_TYPE::GHOST_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
          odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::GHOST_FACE};
      iface.mesh_data.number_of_ghost_cells = 2;
      iface.mesh_data.next_cell_id = {0, 4, 4, 1};
      iface.mesh_data.ghost_cell_global_id = {0, 3};
      iface.mesh_data.ghost_cell_proc = {0, 2};
    } else {
      iface.mesh_data.cell_position = {0.25, 0.75, 0.0, 0.75, 0.75, 0.0};
      iface.mesh_data.cell_size = {0.5, 0.5, 0.0, 0.5, 0.5, 0.0};
      iface.mesh_data.cell_global_id = {2, 3};
      iface.mesh_data.face_types = {
          odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
          odd_solver::FACE_TYPE::GHOST_FACE,    odd_solver::FACE_TYPE::BOUNDARY_FACE,
          odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
          odd_solver::FACE_TYPE::GHOST_FACE,    odd_solver::FACE_TYPE::BOUNDARY_FACE};
      iface.mesh_data.number_of_ghost_cells = 2;
      iface.mesh_data.next_cell_id = {4, 1, 0, 4, 0, 4, 1, 4};
      iface.mesh_data.ghost_cell_global_id = {0, 1};
      iface.mesh_data.ghost_cell_proc = {0, 1};
    }
  }
}

void Test_3D_Interface_Builder(odd_solver::Interface_Data &iface, const bool dd = false) {
  // Define mesh data
  // 8 zones with dx=0.5 dy=0.5 and dz=0.5
  //           front        back
  //          ___  ___    ___  ___
  //         | 2 || 3 |  | 6 || 7 |
  //          ===  ===    ===  ===
  //         | 0 || 1 |  | 4 || 5 |
  //          ---  ---    ---  ---
  // if dd -> cells_node_0={0} cell_node_1={1,2,3] cell_node_2={4,5,6,7}
  size_t n_global_cells = 8;
  size_t n_local_cells = n_global_cells;
  if (dd) {
    Insist(rtt_c4::nodes() == 3, "2D DD test mesh only supports 3 ranks");
    if (rtt_c4::node() == 0)
      n_local_cells = 1;
    else if (rtt_c4::node() == 1)
      n_local_cells = 3;
    else
      n_local_cells = 4;
  }
  iface.mesh_data.domain_decomposed = dd;
  iface.mesh_data.number_of_local_cells = n_local_cells;
  iface.mesh_data.number_of_global_cells = n_global_cells;
  iface.mesh_data.n_dims = 3;
  iface.mesh_data.coord_sys = odd_solver::COORDINATE_SYSTEM::CARTESIAN;
  if (!dd) {
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
  } else {
    if (rtt_c4::node() == 0) {
      iface.mesh_data.cell_position = {0.25, 0.25, 0.25};
      iface.mesh_data.cell_size = {0.5, 0.5, 0.5};
      iface.mesh_data.cell_global_id = {0};
      iface.mesh_data.face_types = {
          odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::GHOST_FACE,
          odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::GHOST_FACE,
          odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::GHOST_FACE};
      iface.mesh_data.number_of_ghost_cells = 3;
      iface.mesh_data.next_cell_id = {1, 0, 1, 1, 1, 2};
      iface.mesh_data.ghost_cell_global_id = {1, 2, 4};
      iface.mesh_data.ghost_cell_proc = {1, 1, 2};
    } else if (rtt_c4::node() == 1) {
      iface.mesh_data.cell_position = {0.75, 0.25, 0.25, 0.25, 0.75, 0.25, 0.75, 0.75, 0.25};
      iface.mesh_data.cell_size = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
      iface.mesh_data.cell_global_id = {1, 2, 3};
      iface.mesh_data.face_types = {
          odd_solver::FACE_TYPE::GHOST_FACE,    odd_solver::FACE_TYPE::BOUNDARY_FACE,
          odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
          odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::GHOST_FACE,
          odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
          odd_solver::FACE_TYPE::GHOST_FACE,    odd_solver::FACE_TYPE::BOUNDARY_FACE,
          odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::GHOST_FACE,
          odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
          odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
          odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::GHOST_FACE};
      iface.mesh_data.number_of_ghost_cells = 4;
      iface.mesh_data.next_cell_id = {0, 3, 3, 2, 3, 1, 3, 2, 0, 3, 3, 2, 1, 3, 0, 3, 3, 3};
      iface.mesh_data.ghost_cell_global_id = {0, 5, 6, 7};
      iface.mesh_data.ghost_cell_proc = {0, 2, 2, 2};
    } else {
      iface.mesh_data.cell_position = {0.25, 0.25, 0.75, 0.75, 0.25, 0.75,
                                       0.25, 0.75, 0.75, 0.75, 0.75, 0.75};
      iface.mesh_data.cell_size = {0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5, 0.5};
      iface.mesh_data.cell_global_id = {4, 5, 6, 7};
      iface.mesh_data.face_types = {
          odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
          odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
          odd_solver::FACE_TYPE::GHOST_FACE,    odd_solver::FACE_TYPE::BOUNDARY_FACE,
          odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
          odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
          odd_solver::FACE_TYPE::GHOST_FACE,    odd_solver::FACE_TYPE::BOUNDARY_FACE,
          odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
          odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
          odd_solver::FACE_TYPE::GHOST_FACE,    odd_solver::FACE_TYPE::BOUNDARY_FACE,
          odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
          odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE,
          odd_solver::FACE_TYPE::GHOST_FACE,    odd_solver::FACE_TYPE::BOUNDARY_FACE};
      iface.mesh_data.number_of_ghost_cells = 4;
      iface.mesh_data.next_cell_id = {4, 1, 4, 2, 0, 4, 0, 4, 4, 3, 1, 4,
                                      4, 3, 0, 4, 2, 4, 2, 4, 1, 4, 3, 4};
      iface.mesh_data.ghost_cell_global_id = {0, 1, 2, 3};
      iface.mesh_data.ghost_cell_proc = {0, 1, 1, 1};
    }
  }
}

} // namespace odd_solver_test

#endif // solver_test_Test_Interface_Builder

//------------------------------------------------------------------------------------------------//
// end of <pkg>/<class>.hh
//------------------------------------------------------------------------------------------------//
