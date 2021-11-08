//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   api/Arguments.cc
 * \author Mathew Cleveland
 * \date   October 26th 2021
 * \brief  Implementation of api arguments
 * \note   Copyright (C) 2021 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "Arguments.hh"
#include "ds++/dbc.hh"

Control_Data::Control_Data() : opacity_file("") {}

void Control_Data::check_arguments() {
  Insist(opacity_file != "", "Opacity file was not specified");
}

Zonal_Data::Zonal_Data()
    : number_of_cells(0), dimensions(0), dx(0.0), dy(0.0), dz(0.0), number_of_mats(0),
      problem_matids(nullptr), number_of_cell_mats(nullptr), cell_mats(nullptr),
      cell_mat_vol_frac(nullptr), cell_mat_density(nullptr) {}

void Zonal_Data::check_arguments() {
  Insist(number_of_cells > 0, "Number of cells must be greater then zero");
  Insist(dimensions > 0, "Dimension must be greater then zero");
  Insist(dx > 0, "dx must be greater then zero");
  Insist(dimensions > 1 ? dy > 0.0 : true, "dy must be greater then zero if dimensions>1");
  Insist(dimensions > 2 ? dz > 0.0 : true, "dz must be greater then zero if dimensions>1");
  Insist(number_of_mats > 0, "Must have more then zero materials");
  Insist(problem_matids != nullptr, "Vector of Problem material IDs was not specified");
  Insist(number_of_cell_mats != nullptr, "Number of materials per cell was not specified");
  Insist(cell_mats != nullptr, "Vector of cell material IDs was not specified");
  Insist(cell_mat_vol_frac != nullptr,
         "Vector of cell material volume fractions was not specified");
  Insist(cell_mat_temperature != nullptr, "Vector of cell material temperatures was not specified");
  Insist(cell_mat_density != nullptr, "Vector of cell material densities was not specified");
}

Output_Data::Output_Data() : ave_opacity_data(nullptr) {}

void Output_Data::check_arguments() {
  Insist(ave_opacity_data != nullptr, "Opacity_data field must be specified");
}

Arguments::Arguments() : control_data(), zonal_data(), output_data() {}

//------------------------------------------------------------------------------------------------//
// end of Arguments.cc
//------------------------------------------------------------------------------------------------//
