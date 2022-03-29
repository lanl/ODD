//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   solver/Interface_Data.hh
 * \author Mathew Cleveland
 * \brief  Define the Solver Interface data class that holds the fundamental interface data in
 * convenient c++ data types (ie vectors rather then pointers)
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef solver_Interface_Data_hh
#define solver_Interface_Data_hh
#include "units/PhysicalConstexprs.hh"
#include <array>
#include <string>
#include <vector>

namespace odd_solver {

//! Face types used to navigate the mesh
enum COORDINATE_SYSTEM { CARTESIAN, CYLINDRICAL, SPHERICAL, N_COORD_TYPES };

//! Face types used to navigate the mesh
enum FACE_TYPE { INTERNAL_FACE, BOUNDARY_FACE, GHOST_FACE, N_FACE_TYPES };

//! Control Data
struct Control_Data {
  bool multigroup{false};
  std::array<bool, 6> reflect_bnd{true, true, true, true, true, true};
  std::array<double, 6> bnd_temp{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
};

//! Mesh Data
struct Mesh_Data {
  bool domain_decomposed{false};
  size_t n_dims{0};
  size_t coord_sys{COORDINATE_SYSTEM::N_COORD_TYPES};
  size_t number_of_local_cells{0};
  size_t number_of_ghost_cells{0};
  size_t number_of_global_cells{0};
  // size=n_cells*3
  std::vector<double> cell_position;
  std::vector<double> cell_size;
  // size=n_cells
  std::vector<size_t> cell_global_id;
  // sieze=n_ghost_cells
  std::vector<size_t> ghost_cell_global_id;
  std::vector<size_t> ghost_cell_proc;
  // size=n_cells*n_faces_per_cell
  std::vector<size_t> face_types;
  // size=n_cells*n_faces_per_cell
  std::vector<size_t> next_cell_id;
};

//! Raw Material data arrays
struct Mat_Data {
  std::string ipcress_filename;
  size_t number_of_mats{0};
  std::vector<size_t> problem_matids;
  std::vector<size_t> number_of_cell_mats;
  std::vector<std::vector<size_t>> cell_mats;
  std::vector<std::vector<double>> cell_mat_vol_frac;
  std::vector<std::vector<double>> cell_mat_temperature;
  std::vector<std::vector<double>> cell_mat_density;
  std::vector<std::vector<double>> cell_mat_specific_heat;
  std::vector<double> cell_rad_eden;
  std::vector<std::array<double, 3>> cell_velocity;
};

//! Raw output data
struct Output_Data {
  std::vector<std::vector<double>> cell_mat_dedv;
  std::vector<double> cell_rad_eden;
};

//! Raw matrix data and the solution vectors
struct Solver_Data {
  // Matrix Data
  std::vector<double> diagonal;
  std::vector<std::vector<double>> off_diagonal;
  std::vector<std::vector<size_t>> off_diagonal_id;
  std::vector<std::vector<size_t>> face_type;
  std::vector<double> source;
  // Extra conservation data;
  std::array<double, 6> boundary_source;
  std::array<double, 6> boundary_leakage;
  // Mat Data
  std::vector<double> cell_density;
  std::vector<double> cell_cve;
  // Initial Conditions
  std::vector<double> cell_eden0;
  std::vector<std::vector<double>> face_flux0;
  std::vector<double> cell_temperature0;
  // Solution Data
  std::vector<double> cell_eden;
  std::vector<std::vector<double>> face_flux;
  std::vector<double> cell_temperature;
  // Ghost Data
  std::vector<double> ghost_cell_temperature;
  std::vector<double> ghost_face_D;
  std::vector<double> ghost_dist_center_to_face;
  std::vector<double> ghost_cell_eden;
};

//! Raw multigroup matrix data and the solution vectors
struct MG_Solver_Data {
  // Matrix Data
  std::vector<std::vector<double>> diagonal;
  std::vector<std::vector<std::vector<double>>> off_diagonal;
  std::vector<size_t> off_diagonal_id;
  std::vector<std::vector<double>> source;
  // Mat Data
  std::vector<double> cell_density;
  std::vector<double> cell_cve;
  // Initial Conditions
  std::vector<std::vector<double>> cell_eden0;
  std::vector<std::vector<std::vector<double>>> face_flux0;
  std::vector<double> cell_temperature0;
  // Solution Data
  std::vector<std::vector<double>> cell_eden;
  std::vector<std::vector<std::vector<double>>> face_flux;
  std::vector<double> cell_temperature;
};

//================================================================================================//
/*!
 * \class Interface_Data
 * \brief
 *
 * This holds the fundamental interface data in simple C++ data types
 *
 */
//================================================================================================//

class Interface_Data {
public:
  //! Default constructors.
  Interface_Data() = default;

  bool valid();

  Control_Data control_data;
  Mesh_Data mesh_data;
  Mat_Data mat_data;
  Output_Data output_data;
};

} // namespace odd_solver

#endif // solver_Interface_Data_hh

//------------------------------------------------------------------------------------------------//
// end of solver/Interface_Data.hh
//------------------------------------------------------------------------------------------------//
