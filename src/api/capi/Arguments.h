//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   api/Arguments.hh
 * \author Mathew Cleveland
 * \brief  C Interface
 * \note   Copyright (C) 2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#ifndef api_arguments_h
#define api_arguments_h

#include <stdio.h>

//! Control data
struct Control_Data {
  char *opacity_file;
  // Constructor to initialize data
  double dt;
  size_t max_iter;
  double min_tol;
  double *bnd_temp;
  double *reflect_bnd;
};

//! Zonal data
struct Zonal_Data {
  // Mesh data
  size_t domain_decomposed;
  size_t dimensions;
  size_t coord_sys;
  size_t number_of_local_cells;
  size_t number_of_ghost_cells;
  size_t number_of_global_cells;
  // cell position for each dimension size=n_cells*3
  double *cell_position;
  // cell size in each dimension size=n_cells*3
  double *cell_size;
  // size=n_cells
  size_t *cell_global_id;
  // cell face type for each face of each cell (0=internal, 1=boundary, 2=ghost)
  // size=n_cells*n_faces_per_cell
  size_t *face_type;
  // local id of the next cell for each face of each cell size=n_cells*n_faces_per_cell
  size_t *next_cell_id;
  // size=n_ghost_cells
  size_t *ghost_cell_global_id;
  size_t *ghost_cell_proc;

  // Material data
  // total number of mats in the problem
  size_t number_of_mats;
  // material id for each mat in the problem
  size_t *problem_matids;
  // number of materials in each cell
  size_t *number_of_cell_mats;
  // list of materials in each cell (strided by cells*mats)
  size_t *cell_mats;
  // void fraction of each material in a cell (strided by cells*mats)
  double *cell_mat_vol_frac;
  // temperature of each material in each cell (strided by cells*mats) [keV]
  double *cell_mat_temperature;
  // density of each material in each cell (strided by cells*mats) [g/cc]
  double *cell_mat_density;
  // cell specific heat [jerks/keV/g]
  double *cell_mat_specific_heat;

  // cell velocity ncells*3 [cm/sh]
  double *cell_velocity;
  // cell energy density [jerks/cc]
  double *cell_erad;
};

//! Output data
struct Output_Data {
  // The cell radiation energy density [jerks/cc]
  double *cell_erad;
  // The cell radiation temperature [keV]
  double *cell_Trad;
  // The change in material energy (strided by cells*mats) [jerks/cc]
  double *cell_mat_delta_e;
};

//! Arguments data
struct Arguments {
  struct Control_Data control_data;
  struct Zonal_Data zonal_data;
  struct Output_Data output_data;
};

typedef struct Arguments Arguments_t;

// Available C functions
void Odd_Diffusion_Solve(Arguments_t *arg);

#endif // api_arguments_h

//------------------------------------------------------------------------------------------------//
// end of api/arguments.hh
//------------------------------------------------------------------------------------------------//
