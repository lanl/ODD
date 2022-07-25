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
#include "odd_release/Release.hh"
#include "solver/Interface_Data.hh"
#include "c4/ParallelUnitTest.hh"
#include "ds++/Release.hh"
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
  std::string file = "not_empty";
  arg.control_data.opacity_file = &file[0];
  arg.control_data.dt = 0.1;
  arg.control_data.max_iter = 1;
  arg.control_data.min_tol = 1.0e-12;
  std::vector<double> bnd_temp(6, 0.0);
  arg.control_data.bnd_temp = &bnd_temp[0];
  std::vector<size_t> reflect_bnd(6, 1);
  arg.control_data.reflect_bnd = &reflect_bnd[0];
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
  // cell 1 (0.1) cell 2 mat_0 (0.1) and mat 1 (0.01)
  std::vector<double> cell_mat_electron_source{0.1, 0.1, 0.01};
  // cell velocity
  std::vector<double> cell_velocity{0.0, 0.0, 0.0, 0.0, 0.0, 0.0};
  // cell radiation energy density
  std::vector<double> cell_rad_eden{32.0, 23.0};
  // cell radiation source
  std::vector<double> cell_rad_source{32.0, 23.0};
  // face flux for each cell (nfaces*ncells)
  std::vector<double> face_flux{0.0, 1.0, 1.0, 0.0};
  arg.zonal_data.cell_mats = &cell_mats[0];
  arg.zonal_data.cell_mat_vol_frac = &cell_mat_vol_frac[0];
  arg.zonal_data.cell_mat_temperature = &cell_mat_temperature[0];
  arg.zonal_data.cell_mat_density = &cell_mat_density[0];
  arg.zonal_data.cell_mat_specific_heat = &cell_mat_specific_heat[0];
  arg.zonal_data.cell_mat_electron_source = &cell_mat_electron_source[0];
  arg.zonal_data.cell_velocity = &cell_velocity[0];
  arg.zonal_data.cell_erad = &cell_rad_eden[0];
  arg.zonal_data.cell_rad_source = &cell_rad_source[0];
  arg.zonal_data.face_flux = &face_flux[0];

  // Check that the pointer matches the data
  FAIL_IF_NOT((*arg.zonal_data.problem_matids) == matids[0]);
  arg.zonal_data.problem_matids++;
  FAIL_IF_NOT((*arg.zonal_data.problem_matids) == matids[1]);
  arg.zonal_data.check_arguments();

  try {
    arg.output_data.check_arguments(arg.zonal_data.number_of_local_cells,
                                    arg.zonal_data.number_of_cell_mats,
                                    2 * arg.zonal_data.dimensions);
  } catch (...) {
    std::cout << "Caught the expected Output_Data::check_arguments assertions" << std::endl;
  }
  std::vector<double> cell_erad = {17000.0, 17000.1};
  std::vector<double> cell_Trad = {18000.0, 18000.1};
  std::vector<double> cell_mat_delta_e = {-3.0, 4.0, -5.0};
  std::vector<double> output_face_flux = {0.0, -1.0, -1.0, 0.0};
  arg.output_data.cell_erad = &cell_erad[0];
  arg.output_data.cell_Trad = &cell_Trad[0];
  arg.output_data.cell_mat_delta_e = &cell_mat_delta_e[0];
  arg.output_data.face_flux = &output_face_flux[0];
  // Check that the pointer matches the data
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(arg.output_data.cell_erad[0], cell_erad[0]));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(arg.output_data.cell_erad[1], cell_erad[1]));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(arg.output_data.cell_Trad[0], cell_Trad[0]));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(arg.output_data.cell_Trad[1], cell_Trad[1]));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(arg.output_data.cell_mat_delta_e[0], cell_mat_delta_e[0]));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(arg.output_data.cell_mat_delta_e[1], cell_mat_delta_e[1]));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(arg.output_data.cell_mat_delta_e[2], cell_mat_delta_e[2]));
  FAIL_IF_NOT(rtt_dsxx::soft_equiv(arg.output_data.face_flux[2], output_face_flux[2]));
  arg.output_data.check_arguments(arg.zonal_data.number_of_local_cells,
                                  arg.zonal_data.number_of_cell_mats,
                                  2 * arg.zonal_data.dimensions);

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "the ODD API arguments seem to be working";
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
void test_dd(rtt_dsxx::UnitTest &ut) {
  FAIL_IF_NOT(rtt_c4::nodes() == 2);
  Arguments arg;
  // Should fail
  try {
    arg.control_data.check_arguments();
  } catch (...) {
    std::cout << "Caught the expected Control_Data::check_arguments assertions" << std::endl;
  }
  // fill control data with "valid junk" and check arguments
  std::string file = "not_empty";
  arg.control_data.opacity_file = &file[0];
  arg.control_data.dt = 0.1;
  arg.control_data.max_iter = 1;
  arg.control_data.min_tol = 1.0e-12;
  std::vector<double> bnd_temp(6, 0.0);
  arg.control_data.bnd_temp = &bnd_temp[0];
  std::vector<size_t> reflect_bnd(6, 1);
  arg.control_data.reflect_bnd = &reflect_bnd[0];
  arg.control_data.check_arguments();

  try {
    arg.zonal_data.check_arguments();
  } catch (...) {
    std::cout << "Caught the expected Zonal_Data::check_arguments assertions" << std::endl;
  }
  // fill the zonal data with "valid junk" and check arguments
  arg.zonal_data.domain_decomposed = 1;
  arg.zonal_data.number_of_local_cells = 1;
  arg.zonal_data.number_of_global_cells = 2;
  arg.zonal_data.number_of_ghost_cells = 1;
  arg.zonal_data.dimensions = 1;
  arg.zonal_data.coord_sys = 0;
  std::vector<double> cell_position;
  std::vector<double> cell_size;
  std::vector<size_t> cell_global_id;
  std::vector<size_t> face_type;
  std::vector<size_t> next_cell_id;
  std::vector<size_t> ghost_cell_global_id;
  std::vector<size_t> ghost_cell_proc;
  // 2 zones | 0 || 1 | with dx=0.5 dy=0 and dz=0
  if (rtt_c4::node() == 0) {
    cell_position = {0.25, 0.0, 0.0};
    cell_size = {0.5, 0.0, 0.0};
    cell_global_id = {0};
    face_type = {odd_solver::FACE_TYPE::BOUNDARY_FACE, odd_solver::FACE_TYPE::GHOST_FACE};
    next_cell_id = {1, 0};
    ghost_cell_global_id = {1};
    ghost_cell_proc = {1};
  }
  if (rtt_c4::node() == 1) {
    cell_position = {0.75, 0.0, 0.0};
    cell_size = {0.5, 0.0, 0.0};
    cell_global_id = {1};
    face_type = {odd_solver::FACE_TYPE::GHOST_FACE, odd_solver::FACE_TYPE::BOUNDARY_FACE};
    next_cell_id = {0, 1};
    ghost_cell_global_id = {0};
    ghost_cell_proc = {0};
  }
  arg.zonal_data.cell_position = &cell_position[0];
  arg.zonal_data.cell_size = &cell_size[0];
  arg.zonal_data.cell_global_id = &cell_global_id[0];
  arg.zonal_data.face_type = &face_type[0];
  arg.zonal_data.next_cell_id = &next_cell_id[0];
  arg.zonal_data.ghost_cell_global_id = &ghost_cell_global_id[0];
  arg.zonal_data.ghost_cell_proc = &ghost_cell_proc[0];

  // global material data
  std::vector<size_t> matids = {19000, 190001};
  arg.zonal_data.number_of_mats = 2;
  arg.zonal_data.problem_matids = &matids[0];
  std::vector<size_t> cell_number_of_mats;
  std::vector<size_t> cell_mats;
  std::vector<double> cell_mat_vol_frac;
  std::vector<double> cell_mat_temperature;
  std::vector<double> cell_mat_density;
  std::vector<double> cell_mat_specific_heat;
  std::vector<double> cell_mat_electron_source;
  std::vector<double> cell_velocity;
  std::vector<double> cell_rad_eden;
  std::vector<double> cell_rad_source;
  std::vector<double> face_flux;
  if (rtt_c4::node() == 0) {
    // cell wise material data
    // 1 material in cell 1 and 2 materials in cell 2
    cell_number_of_mats = {1};
    // cell 1 (mat 1) cell 2 (mat 1 and mat 2)
    cell_mats = {0};
    // cell 1 (1.0) cell 2 mat 1 (0.5) and mat 2 (0.5)
    cell_mat_vol_frac = {1.0};
    // cell 1 (1.0) cell 2 mat 1 (10.0) and mat 2 (1.0)
    cell_mat_temperature = {1.0};
    // cell 1 (10) cell 2 mat 1 (1.0) and mat 2 (10.0)
    cell_mat_density = {10.0};
    // cell 1 (0.1) cell 2 mat 0 (0.1) and mat 1 (0.01)
    cell_mat_specific_heat = {0.1};
    // cell 1 (0.1) cell 2 mat 0 (0.1) and mat 1 (0.01)
    cell_mat_electron_source = {0.1};
    // cell velocity
    cell_velocity = {0.0, 0.0, 0.0};
    // cell radiation energy density
    cell_rad_eden = {32.0};
    // cell radiation source
    cell_rad_source = {32.0};
    // flux on each cell face
    face_flux = {0.0, 1.0};
  }
  if (rtt_c4::node() == 1) {
    // cell wise material data
    // 1 material in cell 1 and 2 materials in cell 2
    cell_number_of_mats = {2};
    // cell 1 (mat 1) cell 2 (mat 1 and mat 2)
    cell_mats = {0, 1};
    // cell 1 (1.0) cell 2 mat 1 (0.5) and mat 2 (0.5)
    cell_mat_vol_frac = {0.5, 0.5};
    // cell 1 (1.0) cell 2 mat 1 (10.0) and mat 2 (1.0)
    cell_mat_temperature = {10.0, 1.0};
    // cell 1 (10) cell 2 mat 1 (1.0) and mat 2 (10.0)
    cell_mat_density = {1.0, 10.0};
    // cell 1 (0.1) cell 2 mat 0 (0.1) and mat 1 (0.01)
    cell_mat_specific_heat = {0.1, 0.01};
    // cell 1 (0.1) cell 2 mat 0 (0.1) and mat 1 (0.01)
    cell_mat_electron_source = {0.1, 0.01};
    // cell velocity
    cell_velocity = {0.0, 0.0, 0.0};
    // cell radiation energy density
    cell_rad_eden = {23.0};
    // cell radiation source
    cell_rad_source = {23.0};
    // flux on each cell face
    face_flux = {1.0, 0.0};
  }
  arg.zonal_data.number_of_cell_mats = &cell_number_of_mats[0];
  arg.zonal_data.cell_mats = &cell_mats[0];
  arg.zonal_data.cell_mat_vol_frac = &cell_mat_vol_frac[0];
  arg.zonal_data.cell_mat_temperature = &cell_mat_temperature[0];
  arg.zonal_data.cell_mat_density = &cell_mat_density[0];
  arg.zonal_data.cell_mat_specific_heat = &cell_mat_specific_heat[0];
  arg.zonal_data.cell_mat_electron_source = &cell_mat_electron_source[0];
  arg.zonal_data.cell_velocity = &cell_velocity[0];
  arg.zonal_data.cell_erad = &cell_rad_eden[0];
  arg.zonal_data.cell_rad_source = &cell_rad_source[0];
  arg.zonal_data.face_flux = &face_flux[0];

  arg.zonal_data.check_arguments();

  try {
    arg.output_data.check_arguments(arg.zonal_data.number_of_local_cells,
                                    arg.zonal_data.number_of_cell_mats,
                                    2 * arg.zonal_data.dimensions);
  } catch (...) {
    std::cout << "Caught the expected Output_Data::check_arguments assertions" << std::endl;
  }
  std::vector<double> cell_erad;
  std::vector<double> cell_Trad;
  std::vector<double> cell_mat_delta_e;
  std::vector<double> output_face_flux;
  if (rtt_c4::node() == 0) {
    cell_erad = {17000.0};
    cell_Trad = {18000.0};
    cell_mat_delta_e = {-3.0};
    output_face_flux = {0.0, -1.0};
    arg.output_data.cell_erad = &cell_erad[0];
    arg.output_data.cell_Trad = &cell_Trad[0];
    arg.output_data.cell_mat_delta_e = &cell_mat_delta_e[0];
    arg.output_data.face_flux = &output_face_flux[0];
  }
  if (rtt_c4::node() == 1) {
    cell_erad = {17000.1};
    cell_Trad = {18000.1};
    cell_mat_delta_e = {4.0, -5.0};
    output_face_flux = {-1.0, 0.0};
    arg.output_data.cell_erad = &cell_erad[0];
    arg.output_data.cell_Trad = &cell_Trad[0];
    arg.output_data.cell_mat_delta_e = &cell_mat_delta_e[0];
    arg.output_data.face_flux = &output_face_flux[0];
  }
  arg.output_data.check_arguments(arg.zonal_data.number_of_local_cells,
                                  arg.zonal_data.number_of_cell_mats,
                                  2 * arg.zonal_data.dimensions);

  if (ut.numFails == 0) {
    std::ostringstream m;
    m << "the ODD API DD arguments seem to be working";
    PASSMSG(m.str());
  }
}

//------------------------------------------------------------------------------------------------//
int main(int argc, char *argv[]) {
  rtt_c4::ParallelUnitTest ut(argc, argv, rtt_odd::release);
  try {
    // >>> UNIT TESTS
    test(ut);
    if (rtt_c4::nodes() == 2) {
      test_dd(ut);
    }
  }
  UT_EPILOG(ut);
}

//------------------------------------------------------------------------------------------------//
// end of tstArgumetns.cc
//------------------------------------------------------------------------------------------------//
