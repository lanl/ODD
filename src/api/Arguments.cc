//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   api/Arguments.cc
 * \author Mathew Cleveland
 * \date   October 26th 2021
 * \brief  Implementation of api arguments
 * \note   Copyright (C) 2018-2020 Triad National Security, LLC.
 *         All rights reserved. */
//------------------------------------------------------------------------------------------------//

#include "Arguments.hh"
#include "ds++/dbc.hh"

namespace odd_api {

Control_Data::Control_Data() : opacity_file("") {}

void Control_Data::check_arguments() {
  Insist(opacity_file != "", "Opacity file was not specified");
}

Zonal_Data::Zonal_Data()
    : number_of_cells(0), dimensions(0), dx(0.0), dy(0.0), dz(0.0), matid(nullptr) {}

void Zonal_Data::check_arguments() {
  Insist(number_of_cells > 0, "Number of cells must be greater then zero");
  Insist(dimensions > 0, "Dimension must be greater then zero");
  Insist(dx > 0, "dx must be greater then zero");
  Insist(dimensions > 1 ? dy > 0.0 : true, "dy must be greater then zero if dimensions>1");
  Insist(dimensions > 2 ? dz > 0.0 : true, "dz must be greater then zero if dimensions>1");
  Insist(matid != nullptr, "Material ID array mulst be specified");
}

Output_Data::Output_Data() : opacity_data(nullptr) {}

void Output_Data::check_arguments() {
  Insist(opacity_data != nullptr, "Opacity_data field must be specified");
}

Arguments::Arguments() : control_data(), zonal_data(), output_data() {}

} // namespace odd_api

//------------------------------------------------------------------------------------------------//
// end of Arguments.cc
//------------------------------------------------------------------------------------------------//
