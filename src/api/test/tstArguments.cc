//--------------------------------------------*-C++-*---------------------------------------------//
/*!
 * \file   api/test/tstArguments.cc
 * \author Mathew Cleveland
 * \date   October 21st 2021
 * \brief  Testing opacity reader
 * \note   Copyright (C) 2021-2022 Triad National Security, LLC., All rights reserved.
 */
//------------------------------------------------------------------------------------------------//

#include "api/Arguments.hh"
#include "solver/Interface_Data.hh"
#include "ds++/Release.hh"
#include "ds++/ScalarUnitTest.hh"
#include "ds++/dbc.hh"

using namespace rtt_dsxx;

//------------------------------------------------------------------------------------------------//
// TESTS
//------------------------------------------------------------------------------------------------//
void test(rtt_dsxx::UnitTest &ut) {
  Arguments arg;
  // Should fail
  try {
    arg.control_data.check_arguments();
  } catch (...) {
    std::cout << "Caught the expected Control_Data::check_arguments assertions" << std::endl;
  }
  // fill control data with "valid junk" and check arguments
  arg.control_data.opacity_file = "not_empty";
  arg.control_data.check_arguments();

  try {
    arg.zonal_data.check_arguments();
  } catch (...) {
    std::cout << "Caught the expected Zonal_Data::check_arguments assertions" << std::endl;
  }
  // fill the zonal data with "valid junk" and check arguments
  arg.zonal_data.domain_decomposed = 0;
  arg.zonal_data.number_of_local_cells = 2;
  arg.zonal_data.number_of_global_cells = 2;
  arg.zonal_data.dimensions = 1;
  arg.zonal_data.coord_sys = 0;
  // 2 zones | 0 || 1 | with dx=0.5 dy=0 and dz=0
  std::vector<double> cell_position{0.25, 0.0, 0.0, 0.75, 0.0, 0.0};
  arg.zonal_data.cell_position = &cell_position[0];
  std::vector<double> cell_size{0.5, 0.0, 0.0, 0.5, 0.0, 0.0};
  arg.zonal_data.cell_size = &cell_size[0];
  std::vector<size_t> cell_global_id{0, 1};
  arg.zonal_data.cell_global_id = &cell_global_id[0];
  std::vector<size_t> face_type{
      odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::INTERNAL_FACE,
      odd_solver::FACE_TYPE::INTERNAL_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE};
  arg.zonal_data.face_type = &face_type[0];
  std::vector<size_t> next_cell_id{2, 1, 0, 2};
  arg.zonal_data.next_cell_id = &next_cell_id[0];

  // global material data
  std::vector<size_t> matids = {19000, 190001};
  arg.zonal_data.number_of_mats = 2;
  arg.zonal_data.problem_matids = &matids[0];
  // cell wise material data
  // 1 material in cell 1 and 2 materials in cell 2
  std::vector<size_t> cell_number_of_mats{1, 2};
  arg.zonal_data.number_of_cell_mats = &cell_number_of_mats[0];
  // cell 1 (mat 1) cell 2 (mat 1 and mat 2)
  std::vector<size_t> cell_mats{0, 0, 1};
  // cell 1 (1.0) cell 2 mat 1 (0.5) and mat 2 (0.5)
  std::vector<double> cell_mat_vol_frac{1.0, 0.5, 0.5};
  // cell 1 (1.0) cell 2 mat 1 (10.0) and mat 2 (1.0)
  std::vector<double> cell_mat_temperature{1.0, 10.0, 1.0};
  // cell 1 (10) cell 2 mat 1 (1.0) and mat 2 (10.0)
  std::vector<double> cell_mat_density{10.0, 1.0, 10.0};
  // cell 1 (0.1) cell 2 mat 0 (0.1) and mat 1 (0.01)
  std::vector<double> cell_mat_specific_heat{0.1, 0.1, 0.01};
  // cell velocity
  std::vector<double> cell_velocity{0.0, 0.0};
  arg.zonal_data.cell_mats = &cell_mats[0];
  arg.zonal_data.cell_mat_vol_frac = &cell_mat_vol_frac[0];
  arg.zonal_data.cell_mat_temperature = &cell_mat_temperature[0];
  arg.zonal_data.cell_mat_density = &cell_mat_density[0];
  arg.zonal_data.cell_mat_specific_heat = &cell_mat_specific_heat[0];
  arg.zonal_data.cell_velocity = &cell_velocity[0];

  // Check that the pointer matches the data
  FAIL_IF_NOT((*arg.zonal_data.problem_matids) == matids[0]);
  arg.zonal_data.problem_matids++;
  FAIL_IF_NOT((*arg.zonal_data.problem_matids) == matids[1]);
  arg.zonal_data.check_arguments();

  try {
    arg.output_data.check_arguments();
  } catch (...) {
    std::cout << "Caught the expected Output_Data::check_arguments assertions" << std::endl;
  }
  std::vector<double> opacity_data = {17000.0, 17000.1};
  arg.output_data.ave_opacity_data = &opacity_data[0];
  // Check that the pointer matches the data
  FAIL_IF_NOT(rtt_dsxx::soft_equiv((*arg.output_data.ave_opacity_data), opacity_data[0]));
  arg.output_data.ave_opacity_data++;
  FAIL_IF_NOT(rtt_dsxx::soft_equiv((*arg.output_data.ave_opacity_data), opacity_data[1]));
  arg.output_data.check_arguments();

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "the ODD API arguments seem to be working";
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  ScalarUnitTest ut(argc, argv, release);
  try {
    // >>> UNIT TESTS
    test(ut);
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstArgumetns.cc
//------------------------------------------------------------------------------------------------//
