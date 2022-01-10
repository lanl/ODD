//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   api/Arguments.cc
 * \author Mathew Cleveland
 * \date   October 26th 2021
 * \brief  Implementation of api arguments
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Arguments.hh"
#include "c4/global.hh"
#include "ds++/DracoMath.hh"
#include "ds++/dbc.hh"

Control_Data::Control_Data() : opacity_file("") {}

//================================================================================================//
/*!
 * \brief Check control data arguments
 *
 */
//================================================================================================//
void Control_Data::check_arguments() {
  Insist(opacity_file != "", "Opacity file was not specified");
}

//================================================================================================//
/*!
 * \brief Check zone data arguments
 *
 */
//================================================================================================//
void Zonal_Data::check_arguments() {
  // Check mesh data first
  Insist(domain_decomposed == 0 || domain_decomposed == 1,
         "Domain decomposed must be true or false (1 or 0)");
  Insist(dimensions > 0 && dimensions < 4, "Dimensions must be (1, 2, or 3)");
  Insist(coord_sys < odd_solver::COORDINATE_SYSTEM::N_COORD_TYPES, "Invalid coordinate system");
  Insist(number_of_local_cells > 0, "Number of local cells must be greater then zero");
  Insist(domain_decomposed == 1 ? number_of_ghost_cells > 0 : true,
         "Number of ghost cells must be greater then zero in domain decomposed problems");
  Insist(number_of_global_cells > 0, "Number of global cells must be greater then zero");
  Insist(domain_decomposed == 0 ? number_of_global_cells == number_of_local_cells : true,
         "Number of global and local cells much match for domain replicated problems");
  Insist(cell_position != nullptr, "Must define the cell position array");
  Insist(cell_size != nullptr, "Must define the cell size array array");
  Insist(cell_global_id != nullptr, "Must define the cell global id array");
  Insist(face_type != nullptr, "Must define face type array");
  Insist(next_cell_id != nullptr, "Must define next cell array");
  // Check local cell array values
  for (size_t i = 0; i < number_of_local_cells; i++) {
    Insist(cell_global_id[i] < number_of_global_cells, "Cell global id must be greater then zero");
    for (size_t d = 0; d < dimensions; d++) {
      Insist(rtt_dsxx::isFinite(cell_position[i * 3 + d]), "Cell position must be finite");
      Insist(rtt_dsxx::isFinite(cell_size[i * 3 + d]), "Cell size must be finite");
      Insist(face_type[i * dimensions * 2 + d * 2 + d] < odd_solver::FACE_TYPE::N_FACE_TYPES,
             "Cell lower face type must be a valid type (0=internal, 1=boundary, or 2=ghost");
      Insist(face_type[i * dimensions * 2 + d * 2 + d + 1] < odd_solver::FACE_TYPE::N_FACE_TYPES,
             "Cell upper face type must be a valid type (0=internal, 1=boundary, or 2=ghost");
      Insist(face_type[i * dimensions * 2 + d * 2 + d] < odd_solver::FACE_TYPE::N_FACE_TYPES,
             "Cell lower face type must be a valid type (0=internal, 1=boundary, or 2=ghost");
      // cell_id == number_of_local_cells is a boundary designator
      Insist(next_cell_id[i * dimensions * 2 + d * 2 + d] <=
                 std::max(number_of_local_cells, number_of_ghost_cells),
             "next lower cell must be valid");
      Insist(next_cell_id[i * dimensions * 2 + d * 2 + d + 1] <=
                 std::max(number_of_local_cells, number_of_ghost_cells),
             "next upper cell must be valid");
    }
  }

  Insist(domain_decomposed == 1 ? ghost_cell_global_id != nullptr : true,
         "Must define the ghost cell global id array for domain decomposed problems");
  Insist(domain_decomposed == 1 ? ghost_cell_proc != nullptr : true,
         "Must define the ghost cell proc array for domain decomposed problems");
  // Check ghost array values
  if (domain_decomposed == 1) {
    for (size_t i = 0; i < number_of_global_cells; i++) {
      Insist(ghost_cell_global_id[i] < number_of_global_cells, "Cell global id must be valid");
      Insist(ghost_cell_proc[i] < static_cast<size_t>(rtt_c4::nodes()),
             "Cell proc must be bound by mpi range");
    }
  }

  // Check mat data second
  Insist(number_of_mats > 0, "Must have more then zero materials");
  Insist(problem_matids != nullptr, "Vector of Problem material IDs was not specified");
  Insist(number_of_cell_mats != nullptr, "Number of materials per cell was not specified");
  Insist(cell_mats != nullptr, "Vector of cell material IDs was not specified");
  Insist(cell_mat_vol_frac != nullptr,
         "Vector of cell material volume fractions was not specified");
  Insist(cell_mat_temperature != nullptr, "Vector of cell material temperatures was not specified");
  Insist(cell_mat_density != nullptr, "Vector of cell material densities was not specified");
  Insist(cell_mat_specific_heat != nullptr,
         "Vector of cell material specific heat was not specified");
  Insist(cell_velocity != nullptr, "Cell velocity was not specified");
  // Check mat data array values
  size_t index = 0;
  for (size_t i = 0; i < number_of_local_cells; i++) {
    Insist(number_of_cell_mats[i] <= number_of_mats,
           "Number of cell mats must be less then or equal to the total number of mats");
    Insist(rtt_dsxx::isFinite(cell_velocity[i]), "The cell velocity must be finite");
    for (size_t m = 0; m < number_of_cell_mats[i]; m++) {
      Insist(cell_mats[index] < number_of_mats,
             "Cell mat must be less then or equal to the total number of mats");
      Insist(!(cell_mat_vol_frac[index] < 0.0) && cell_mat_vol_frac[index] < (1.0 + 1.0e-6),
             "Cell mat vol frac must be bound between zero and one");
      Insist(rtt_dsxx::isFinite(cell_mat_temperature[index]) &&
                 !(cell_mat_temperature[index] < 0.0),
             "Cell mat temperature must be non-negative and finite");
      Insist(rtt_dsxx::isFinite(cell_mat_density[index]) && !(cell_mat_density[index] < 0.0),
             "Cell mat density must be non-negative and finite");
      Insist(rtt_dsxx::isFinite(cell_mat_specific_heat[index]) &&
                 !(cell_mat_specific_heat[index] < 0.0),
             "Cell mat temperature must be non-negative and finite");
      index++;
    }
  }
}

//================================================================================================//
/*!
 * \brief Check output data arguments
 *
 */
//================================================================================================//
void Output_Data::check_arguments() {
  Insist(ave_opacity_data != nullptr, "Opacity_data field must be specified");
}

//================================================================================================//
/*!
 * \brief Default arguments constructor
 *
 */
//================================================================================================//
Arguments::Arguments() : control_data(), zonal_data(), output_data() {}

//------------------------------------------------------------------------------------------------//
// end of Arguments.cc
//------------------------------------------------------------------------------------------------//
