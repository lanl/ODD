//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   api/Arguments.hh
 * \author Mathew Cleveland
 * \brief  API argument definitions 
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef api_arguments_hh
#define api_arguments_hh

#include "solver/Interface_Data.hh"
#include <string>

extern "C" {

//! Control data
struct Control_Data {
  char *opacity_file{nullptr};
  // Constructor to initialize data
  double dt{0.0};
  size_t max_iter{0};
  double min_tol{0};
  Control_Data() = default;
  void check_arguments() const;
};

//! Zonal data
struct Zonal_Data {
  // Mesh data
  size_t domain_decomposed{0};
  size_t dimensions{0};
  size_t coord_sys{odd_solver::COORDINATE_SYSTEM::CARTESIAN};
  size_t number_of_local_cells{0};
  size_t number_of_ghost_cells{0};
  size_t number_of_global_cells{0};
  // cell position for each dimension size=n_cells*3
  double *cell_position{nullptr};
  // cell size in each dimension size=n_cells*3
  double *cell_size{nullptr};
  // size=n_cells
  size_t *cell_global_id{nullptr};
  // cell face type for each face of each cell (0=internal, 1=boundary, 2=ghost)
  // size=n_cells*n_faces_per_cell
  size_t *face_type{nullptr};
  // local id of the next cell for each face of each cell size=n_cells*n_faces_per_cell
  size_t *next_cell_id{nullptr};
  // size=n_ghost_cells
  size_t *ghost_cell_global_id{nullptr};
  size_t *ghost_cell_proc{nullptr};

  // Material data
  // total number of mats in the problem
  size_t number_of_mats{0};
  // material id for each mat in the problem
  size_t *problem_matids{nullptr};
  // number of materials in each cell
  size_t *number_of_cell_mats{nullptr};
  // list of materials in each cell (strided by cells*mats)
  size_t *cell_mats{nullptr};
  // void fraction of each material in a cell (strided by cells*mats)
  double *cell_mat_vol_frac{nullptr};
  // temperature of each material in each cell (strided by cells*mats) [keV]
  double *cell_mat_temperature{nullptr};
  // density of each material in each cell (strided by cells*mats) [g/cc]
  double *cell_mat_density{nullptr};
  // cell specific heat [jerks/keV/g]
  double *cell_mat_specific_heat{nullptr};

  // cell velocity ncells*3 [cm/sh]
  double *cell_velocity{nullptr};
  // cell energy density [jerks/cc]
  double *cell_erad{nullptr};

  // constructor to initialize data
  Zonal_Data() = default;
  void check_arguments() const;
};

//! Output data
struct Output_Data {

  // The cell radiation energy density [jerks/cc]
  double *cell_erad{nullptr};
  // The cell radiation temperature [keV]
  double *cell_Trad{nullptr};
  // The change in material energy (strided by cells*mats) [jerks/cc]
  double *cell_mat_delta_e{nullptr};

  // constructor to initialize data
  Output_Data() = default;
  void check_arguments(const size_t ncells, const size_t *number_of_cell_mats) const;
};

//! Arguments data
struct Arguments {
  Control_Data control_data;
  Zonal_Data zonal_data;
  Output_Data output_data;
  // constructor to initialize all data fields
  Arguments();
};

} // extern "C"

#endif // api_arguments_hh

//------------------------------------------------------------------------------------------------//
// end of api/arguments.hh
//------------------------------------------------------------------------------------------------//
